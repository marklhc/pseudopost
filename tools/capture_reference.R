#!/usr/bin/env Rscript
# =============================================================================
# tools/capture_reference.R
#
# Reference capture: run the ORIGINAL pairwise-cfa reference implementation
# (pairwise_func.R + sim_shared_utils.R + my_sim-1fcfa.R wiring) on a small
# fixed ordinal dataset, and record est/se/q5/q95 per parameter per method in
#
#     tests/testthat/_reference/expectations.csv
#
# so the refactored `pseudopost` package can be regression-tested against the
# reference behavior. This script is a development tool; it is NOT part of the
# package (see .Rbuildignore) and does not modify anything under
# /home/marklai/pairwise-cfa.
#
# Methods captured (labels): NoAdj, IJSE, IJSE2, CurvAdj, CurvAdj2, CurvAdj3,
# MagAdj.
# =============================================================================

suppressMessages({
  library(cmdstanr)
  library(numDeriv)
  library(posterior)
})

REF_DIR <- "/home/marklai/pairwise-cfa"
PKG_DIR <- "/home/marklai/pseudopost"
source(file.path(REF_DIR, "sim_shared_utils.R"))
source(file.path(REF_DIR, "pairwise_func.R"))

T_START <- Sys.time()
log_t <- function(...) cat(sprintf("[%6.1f s] ",
                                   as.numeric(Sys.time() - T_START, units = "secs")),
                           ..., "\n")

# -----------------------------------------------------------------------------
# Data: first 40 rows of HolzingerSwineford1939 (x1 math, x2 reading,
# x3 science), each cut into 3 ordered categories via `cut(..., 3)` --
# deterministic EQUAL-WIDTH bins across each column's range (not quantile /
# equal-frequency breaks), no RNG involved. p = 3, k = 3, N = 40.
# -----------------------------------------------------------------------------
hs <- lavaan::HolzingerSwineford1939
n_rows <- 40L
dat_ord <- data.frame(
  x1 = as.ordered(cut(hs$x1[seq_len(n_rows)], 3, include.lowest = TRUE)),
  x2 = as.ordered(cut(hs$x2[seq_len(n_rows)], 3, include.lowest = TRUE)),
  x3 = as.ordered(cut(hs$x3[seq_len(n_rows)], 3, include.lowest = TRUE))
)
p <- 3L
k <- 3L
N <- nrow(dat_ord)
stopifnot(N == 40L, p == 3L, k == 3L)
cat("Pair tables:\n")
print(table(dat_ord$x1, dat_ord$x2))
print(table(dat_ord$x2, dat_ord$x3))

make_param_names <- function(p, k, delta = FALSE) {
  suff <- if (delta) "_delta" else ""
  lam_nms <- paste0("lambda", suff, "[", seq_len(p), "]")
  thr_nms <- unlist(lapply(seq_len(p), function(i)
    paste0("thres", suff, "[", i, ",", seq_len(k - 1), "]")))
  c(lam_nms, thr_nms)
}

# same structure as my_sim-1fcfa.R make_stan_dat
make_stan_dat <- function(p, k, N, dat) {
  pc <- prep_pairwise_counts(dat)
  num_pairs <- if (!is.null(pc$num_pairs)) pc$num_pairs else dim(pc$ns)[1]
  list(p = p, k = k, ns = pc$ns, num_pairs = as.integer(num_pairs), N = N,
       y = sapply(dat, as.integer))
}
stan_dat <- make_stan_dat(p, k, N, dat_ord)

# Reasonable delta-scale starting values so the optimizers converge.
lam0 <- c(0.75, 0.65, 0.55)
tau0 <- list(c(-0.4, 0.6), c(-0.6, 0.4), c(-0.3, 0.5))
stan_init <- function() {
  list(
    lambda_1 = lam0[1] / sqrt(1 - lam0[1]^2),
    lambda_rest = lam0[-1] / sqrt(1 - lam0[-1]^2),
    thres = mapply(function(x, y) x / sqrt(1 - y^2), tau0, lam0,
                   SIMPLIFY = FALSE)
  )
}

pnmsd <- make_param_names(p, k, delta = TRUE)
pnms <- make_param_names(p, k, delta = FALSE)

# -----------------------------------------------------------------------------
# Compile the tinypair fixture (identical to the reference model).
# -----------------------------------------------------------------------------
stan_file <- file.path(PKG_DIR, "inst", "test-models", "tinypair.stan")
log_t("compiling", stan_file)
pcfa_mod <- cmdstan_model(stan_file, force_recompile = TRUE,
                          compile_standalone = TRUE)
log_t("compiled")

# -----------------------------------------------------------------------------
# Main fit (fixed seed).
# -----------------------------------------------------------------------------
mcmc_args <- list(chains = 3L, iter_warmup = 1000L, iter_sampling = 1000L,
                  seed = 42L, refresh = 0, show_messages = FALSE,
                  show_exceptions = FALSE)
log_t("fitting main model (3 x 2000 iters)")
pcfa_fit <- do.call(
  pcfa_mod$sample,
  c(list(data = c(stan_dat, list(magnitude_adj = 1.0, use_priors = 1L)),
         init = stan_init), mcmc_args)
)
log_t("fit done")

# long-form row builder
mk_rows <- function(method, param, est, se, q5, q95) {
  data.frame(method = method, param = param, est = est, se = se, q5 = q5,
             q95 = q95, stringsAsFactors = FALSE)
}
# summarise_draws (posterior) -> 4 stats aligned to pnmsd
rows_from_summary <- function(method, summ_df) {
  idx <- match(pnmsd, summ_df$variable)
  stopifnot(all(!is.na(idx)))
  mk_rows(method, pnmsd, summ_df$mean[idx], summ_df$sd[idx],
          summ_df$q5[idx], summ_df$q95[idx])
}

all_rows <- list()

# ---- NoAdj -----------------------------------------------------------------
summ_noadj <- as.data.frame(pcfa_fit$summary(pnmsd))
all_rows$NoAdj <- rows_from_summary("NoAdj", summ_noadj)
log_t("NoAdj done")

# ---- IJSE (Wald SE from constrained-space IJ covariance) -------------------
# NOTE: positional alignment, exactly as my_sim-1fcfa.R: cov_ij_direct(fit,
# pnmsd) yields the 9 SEs in pnmsd order, and the NoAdj means are in the same
# order. (pnms and pnmsd are name-mismatches; the reference indexes by
# position, not by name.)
pcfa_fit$init_model_methods()
cov_ij1 <- cov_ij_direct(pcfa_fit, pnmsd)
ijse_se <- sqrt(diag(cov_ij1))
ijse_mean <- summ_noadj$mean[match(pnmsd, summ_noadj$variable)]
q95_1sided <- qnorm(0.95)
all_rows$IJSE <- mk_rows(
  "IJSE", pnmsd,
  est = ijse_mean,
  se = ijse_se,
  q5 = ijse_mean - q95_1sided * ijse_se,
  q95 = ijse_mean + q95_1sided * ijse_se
)
log_t("IJSE done")

# ---- IJSE2 (full-draws transform, IJ) --------------------------------------
ij_draws <- transform_draws_ij(pcfa_fit, c("lambda_delta", "thres_delta"))
all_rows$IJSE2 <- rows_from_summary("IJSE2", as.data.frame(posterior::summarise_draws(ij_draws)))
log_t("IJSE2 done")

# ---- CurvAdj (point Godambe at MLE, jacobian = FALSE) ----------------------
opt_mle <- pcfa_mod$optimize(
  data = c(stan_dat, list(magnitude_adj = 1.0, use_priors = 0L)),
  refresh = 0, show_messages = FALSE, init = stan_init
)
log_t("MLE optimize done")
m_godambe_1 <- get_godambe_m(opt_mle)
c_draws_1 <- transform_draws_godambe(pcfa_fit, c("lambda_delta", "thres_delta"), m_godambe_1)
all_rows$CurvAdj <- rows_from_summary("CurvAdj", as.data.frame(posterior::summarise_draws(c_draws_1)))
log_t("CurvAdj done")

# ---- CurvAdj2 (point Godambe at MAP, jacobian = TRUE) ----------------------
opt_map <- pcfa_mod$optimize(
  data = c(stan_dat, list(magnitude_adj = 1.0, use_priors = 1L)),
  refresh = 0, show_messages = FALSE, init = stan_init, jacobian = TRUE
)
log_t("MAP optimize done")
m_godambe_2 <- get_godambe_m(opt_map)
c_draws_2 <- transform_draws_godambe(pcfa_fit, c("lambda_delta", "thres_delta"), m_godambe_2)
all_rows$CurvAdj2 <- rows_from_summary("CurvAdj2", as.data.frame(posterior::summarise_draws(c_draws_2)))
log_t("CurvAdj2 done")

# ---- CurvAdj3 (Godambe with empirical postvar bread) -----------------------
m_godambe_3 <- get_godambe_m(pcfa_fit, htype = "postvar")
c_draws_3 <- transform_draws_godambe(pcfa_fit, c("lambda_delta", "thres_delta"), m_godambe_3)
all_rows$CurvAdj3 <- rows_from_summary("CurvAdj3", as.data.frame(posterior::summarise_draws(c_draws_3)))
log_t("CurvAdj3 done")

# ---- MagAdj (scalar magnitude refit, NoAdj summary of the refit) -----------
mag <- get_mag_adj(opt_mle)
cat("magnitude_adj =", round(mag, 6), "\n")
pcfa_adj_fit <- do.call(
  pcfa_mod$sample,
  c(list(data = c(stan_dat, list(magnitude_adj = mag, use_priors = 1L)),
         init = stan_init), mcmc_args)
)
log_t("MagAdj refit done")
summ_adj <- as.data.frame(pcfa_adj_fit$summary(pnmsd))
all_rows$MagAdj <- rows_from_summary("MagAdj", summ_adj)

# -----------------------------------------------------------------------------
# Write expectations.csv (long form: method, param, est, se, q5, q95).
# -----------------------------------------------------------------------------
expectations <- do.call(rbind, all_rows)
rownames(expectations) <- NULL
out_dir <- file.path(PKG_DIR, "tests", "testthat", "_reference")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, "expectations.csv")
write.csv(expectations, out_file, row.names = FALSE)
log_t("wrote", out_file, "(", nrow(expectations), "rows )")

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat("\n=== capture summary ===\n")
for (m in names(all_rows)) {
  df <- all_rows[[m]]
  vals <- c(df$est, df$se, df$q5, df$q95)
  cat(sprintf("%-9s rows=%2d  finite=%2d/%2d  NA/err=%d\n",
              m, nrow(df), sum(is.finite(vals)), length(vals),
              sum(!is.finite(vals))))
}
cat("\ncapture wall time:",
    round(as.numeric(difftime(Sys.time(), T_START, units = "mins")), 1), "min\n")
