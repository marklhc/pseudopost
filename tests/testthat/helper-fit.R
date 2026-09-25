# ---------------------------------------------------------------------------
# helper-fit.R
#
# Session-level cache for the Stan-backed tests.
#
# COMPILER GOTCHA (why this exists): cmdstanr model methods (`hessian`,
# `constrain_variables`, `unconstrain_draws`) FAIL with "Model methods cannot
# be used with a pre-compiled Stan executable" when the executable is reused
# from the cmdstanr cache (force_recompile = FALSE). They DO work on a fresh
# compile within the same R session. Therefore each fixture is compiled
# exactly ONCE per session with force_recompile = TRUE, and the in-memory
# CmdStanModel object (NOT the exe path) is cached in `.pp_cache` so every
# test file reuses the same live model.
#
# All accessors are lazy: pure-R test files never trigger a compile.
# ---------------------------------------------------------------------------

.pp_cache <- new.env(parent = emptyenv())

#' Skip the Stan-backed integration tests when they cannot run.
#'
#' Skips when (a) PP_SKIP_STAN=1 (fast pure-R iteration), or (b) the CmdStan
#' toolchain the model methods require is not fully available (cmdstanr, Rcpp,
#' and a CmdStan installation). This keeps the package checkable in minimal
#' environments while still exercising the integration tests where the
#' toolchain is present.
skip_stan <- function() {
  if (identical(Sys.getenv("PP_SKIP_STAN"), "1")) {
    skip("stan tests skipped (PP_SKIP_STAN=1)")
  }
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    skip("cmdstanr not installed")
  }
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    skip("Rcpp not installed (required by cmdstanr model methods)")
  }
  if (is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL))) {
    skip("CmdStan toolchain not found (cmdstanr::cmdstan_path())")
  }
}

# Resolve inst/test-models/ (works under test_local / load_all; fall back to
# the source tree for standalone sourcing).
pp_models_dir <- function() {
  d <- system.file("test-models", package = "pseudopost")
  if (!nzchar(d) || !dir.exists(d)) d <- "inst/test-models"
  d
}

# --- truelik: TRUE per-observation likelihood (invariant tests) ------------

get_truelik_model <- function() {
  if (is.null(pp <- .pp_cache$truelik_model)) {
    .pp_cache$truelik_model <- pp <- cmdstanr::cmdstan_model(
      file.path(pp_models_dir(), "truelik.stan"),
      force_recompile = TRUE, quiet = TRUE)
  }
  pp
}

# Fixed pseudo-observations for the truelik fit (seed 42, N = 100).
pp_truelik_y <- function() {
  if (is.null(y <- .pp_cache$truelik_y)) {
    set.seed(42)
    .pp_cache$truelik_y <- y <- rnorm(100, mean = 1.2, sd = 1)
  }
  y
}

pp_truelik_dat <- function() {
  list(N = 100L, y = pp_truelik_y(), magnitude_adj = 1.0, use_priors = 1L)
}

get_truelik_fit <- function() {
  if (is.null(fit <- .pp_cache$truelik_fit)) {
    .pp_cache$truelik_fit <- fit <- get_truelik_model()$sample(
      data = pp_truelik_dat(),
      chains = 2L, iter_warmup = 1000L, iter_sampling = 1000L,
      seed = 42L, refresh = 0, show_messages = FALSE,
      show_exceptions = FALSE)
  }
  fit
}

# --- tinypair: pairwise polychoric CFA (regression vs reference) -----------

# Self-contained replica of the reference `prep_pairwise_counts` +
# `make_stan_dat` (tools/capture_reference.R): first 40 rows of
# HolzingerSwineford1939, each column cut into 3 ordered categories via
# `cut(..., 3)` -- deterministic EQUAL-WIDTH bins across each column's range
# (not quantile/equal-frequency breaks), no RNG. Pair order (1,2),(1,3),(2,3);
# ns[s, r, c] is the row-r / col-c count of pair s, zero-padded to k x k.
pp_make_tinypair_stan_dat <- function() {
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

  pairs <- combn(p, 2L) # columns: (1,2), (1,3), (2,3)
  ns <- array(0L, dim = c(ncol(pairs), k, k))
  for (s in seq_len(ncol(pairs))) {
    i <- pairs[1, s]
    j <- pairs[2, s]
    tbl <- table(dat_ord[[i]], dat_ord[[j]])
    ns[s, seq_len(nrow(tbl)), seq_len(ncol(tbl))] <- as.integer(tbl)
  }
  list(p = p, k = k, ns = ns, num_pairs = as.integer(ncol(pairs)),
       N = N, y = sapply(dat_ord, as.integer))
}

# Reasonable delta-scale starting values (same as tools/capture_reference.R).
pp_tinypair_init <- function() {
  lam0 <- c(0.75, 0.65, 0.55)
  tau0 <- list(c(-0.4, 0.6), c(-0.6, 0.4), c(-0.3, 0.5))
  function() list(
    lambda_1 = lam0[1] / sqrt(1 - lam0[1]^2),
    lambda_rest = lam0[-1] / sqrt(1 - lam0[-1]^2),
    thres = mapply(function(x, y) x / sqrt(1 - y^2), tau0, lam0,
                   SIMPLIFY = FALSE)
  )
}

get_tinypair_model <- function() {
  if (is.null(pp <- .pp_cache$tinypair_model)) {
    .pp_cache$tinypair_model <- pp <- cmdstanr::cmdstan_model(
      file.path(pp_models_dir(), "tinypair.stan"),
      force_recompile = TRUE, quiet = TRUE)
  }
  pp
}

get_tinypair_fit <- function() {
  if (is.null(fit <- .pp_cache$tinypair_fit)) {
    .pp_cache$tinypair_fit <- fit <- get_tinypair_model()$sample(
      data = c(pp_make_tinypair_stan_dat(),
               list(magnitude_adj = 1.0, use_priors = 1L)),
      init = pp_tinypair_init(),
      chains = 3L, iter_warmup = 1000L, iter_sampling = 1000L,
      seed = 42L, refresh = 0, show_messages = FALSE,
      show_exceptions = FALSE)
  }
  fit
}

# Reference parameter labels (9 transformed scalars, capture_reference.R
# `make_param_names(p, k, delta = TRUE)` order).
pp_tinypair_param_names <- function() {
  c("lambda_delta[1]", "lambda_delta[2]", "lambda_delta[3]",
    "thres_delta[1,1]", "thres_delta[1,2]",
    "thres_delta[2,1]", "thres_delta[2,2]",
    "thres_delta[3,1]", "thres_delta[3,2]")
}

# --- nomag: per-observation likelihood WITHOUT a magnitude_adj argument -----
# Exercises the "model cannot be rescaled" path of the magnitude entry points.

get_nomag_model <- function() {
  if (is.null(pp <- .pp_cache$nomag_model)) {
    .pp_cache$nomag_model <- pp <- cmdstanr::cmdstan_model(
      file.path(pp_models_dir(), "nomag.stan"),
      force_recompile = TRUE, quiet = TRUE)
  }
  pp
}

# Fixed pseudo-observations for the nomag fit (seed 7, N = 50).
pp_nomag_y <- function() {
  if (is.null(y <- .pp_cache$nomag_y)) {
    set.seed(7)
    .pp_cache$nomag_y <- y <- rnorm(50, mean = 1.0, sd = 1)
  }
  y
}

pp_nomag_dat <- function() {
  list(N = 50L, y = pp_nomag_y())
}

get_nomag_fit <- function() {
  if (is.null(fit <- .pp_cache$nomag_fit)) {
    .pp_cache$nomag_fit <- fit <- get_nomag_model()$sample(
      data = pp_nomag_dat(),
      chains = 2L, iter_warmup = 200L, iter_sampling = 200L,
      seed = 7L, refresh = 0, show_messages = FALSE,
      show_exceptions = FALSE)
  }
  fit
}
