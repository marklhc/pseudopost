# ---------------------------------------------------------------------------
# Regression vs the reference implementation (tools/capture_reference.R).
# Same machine + seed + model source => reproducible draws, so the package
# table must match tests/testthat/_reference/expectations.csv.
#
# Parameter naming: the comparisons below request the 9 bracketed scalar
# names (pnmsd) directly, exercising constrain_draws()' per-scalar extraction
# for bracketed names (the former bracketed-name bug in R/draws.R is now
# fixed; the last test in this file is the regression guard). fit$draws() reports
# these in the reference pnmsd order, so the package table's `param` column
# matches the reference CSV's `param` column one-to-one.
# ---------------------------------------------------------------------------

pp_ref_label <- c(noadj = "NoAdj", ijse = "IJSE", ijse2 = "IJSE2",
                  curvadj = "CurvAdj", curvadj2 = "CurvAdj2",
                  curvadj3 = "CurvAdj3")

# Compare one column of the package table against the reference rows, with an
# informative worst-case failure message.
pp_check_col <- function(a, r, i, a_col, r_col, tol, label) {
  d <- abs(a[[a_col]] - r[[r_col]][i])
  j <- which.max(d)
  expect_true(max(d) < tol,
              info = sprintf("%s: worst |pkg - ref| = %.6g (tol %g) at param '%s' (pkg %.6g vs ref %.6g)",
                             label, d[j], tol, a$param[j],
                             a[[a_col]][j], r[[r_col]][i][j]))
}

test_that("regression: package matches the reference implementation", {
  skip_stan()
  fit <- get_tinypair_fit()
  mod <- get_tinypair_model()
  stan_dat <- pp_make_tinypair_stan_dat()
  init <- pp_tinypair_init()
  pnmsd <- pp_tinypair_param_names()

  pp <- adjust_pseudo_posterior(
    fit, pnmsd, # bracketed scalar names (see header)
    methods = c("noadj", "ijse", "ijse2", "curvadj", "curvadj2", "curvadj3"),
    model = mod, data = stan_dat, init = init, ci_level = 0.90,
    verbose = FALSE)
  expect_setequal(pp$table$param, pnmsd)
  expect_identical(pp$meta$warnings, character(0))

  ref <- read.csv(test_path("_reference", "expectations.csv"))
  expect_setequal(ref$method, c("NoAdj", "IJSE", "IJSE2", "CurvAdj",
                                "CurvAdj2", "CurvAdj3", "MagAdj"))

  # (curvadj/curvadj2 verified in test-correctness.R: the transform now maps the empirical draw covariance to the target, not the point's inverse-Hessian)
  for (m in c("noadj", "ijse2", "curvadj3")) {
    a <- pp$table[pp$table$method == m, ]
    r <- ref[ref$method == pp_ref_label[[m]], ]
    i <- match(a$param, r$param)
    expect_true(all(!is.na(i)),
                info = paste("params missing from reference:",
                             paste(setdiff(a$param, r$param), collapse = ", ")))
    pp_check_col(a, r, i, "est", "est", 1e-4, m)
    pp_check_col(a, r, i, "se", "se", 5e-3, m)
    pp_check_col(a, r, i, "lo", "q5", 1e-2, m)
    pp_check_col(a, r, i, "hi", "q95", 1e-2, m)
  }

  # ijse: by-design estimator difference. The reference IJSE reports the
  # UNADJUSTED posterior mean as est with a Wald SE from the direct
  # constrained-space covariance; the package ijse reports the mean/SD of the
  # transformed draws (identical transform to ijse2). Hence: se vs IJSE
  # loosely (0.05), est vs the IJSE2 row tightly (same draw transform).
  a <- pp$table[pp$table$method == "ijse", ]
  r_ijse <- ref[ref$method == "IJSE", ]
  r_ijse2 <- ref[ref$method == "IJSE2", ]
  i <- match(a$param, r_ijse$param)
  expect_true(all(!is.na(i)))
  pp_check_col(a, r_ijse, i, "se", "se", 0.05, "ijse vs IJSE (by design)")
  j <- match(a$param, r_ijse2$param)
  expect_true(all(!is.na(j)))
  pp_check_col(a, r_ijse2, j, "est", "est", 1e-4, "ijse est vs IJSE2 est")
})

test_that("regression: MagAdj refit matches the reference", {
  skip_stan()
  fit <- get_tinypair_fit()
  mod <- get_tinypair_model()
  stan_dat <- pp_make_tinypair_stan_dat()
  init <- pp_tinypair_init()
  pnmsd <- pp_tinypair_param_names()

  ma <- refit_magnitude_adjust(
    fit, model = mod, data = stan_dat, init = init,
    params = pnmsd, # bracketed scalar names (see header)
    mcmc_args = list(chains = 3L, iter_warmup = 1000L, iter_sampling = 1000L,
                     seed = 42L, refresh = 0, show_messages = FALSE))
  expect_s3_class(ma, "pseudo_post")
  expect_true(is.finite(ma$magnitude) && ma$magnitude > 0)
  # The refit is summarized as an unadjusted posterior.
  expect_setequal(unique(ma$table$method), "noadj")

  ref <- read.csv(test_path("_reference", "expectations.csv"))
  r <- ref[ref$method == "MagAdj", ]
  a <- ma$table
  i <- match(a$param, r$param)
  expect_true(all(!is.na(i)),
              info = paste("params missing from reference:",
                           paste(setdiff(a$param, r$param), collapse = ", ")))
  pp_check_col(a, r, i, "est", "est", 1e-4, "MagAdj")
  pp_check_col(a, r, i, "se", "se", 5e-3, "MagAdj")
})

test_that("bracketed scalar params are supported by the draw transform", {
  skip_stan()
  fit <- get_tinypair_fit()
  u <- fit$unconstrain_draws(format = "draws_matrix")
  c_naive <- cov(u)
  # T = the draw covariance (C_emp) makes the transform the identity,
  # isolating parameter naming from the covariance mapping.
  sub <- c("lambda_delta[1]", "thres_delta[1,1]")
  # Regression guard for the former bracketed-name bug in constrain_draws()
  # (R/draws.R): bracketed scalar names must run WITHOUT error.
  res_bracketed <- NULL
  expect_no_error(res_bracketed <- transform_to_target(fit, sub,
                                                       T = c_naive))
  expect_true(is.matrix(res_bracketed) && is.numeric(res_bracketed))
  expect_identical(colnames(res_bracketed), sub)
  expect_equal(dim(res_bracketed), c(3000L, 2L))
  # The bracketed result must equal the corresponding columns of the
  # top-level request ...
  res_toplevel <- transform_to_target(fit, c("lambda_delta", "thres_delta"),
                                      T = c_naive)
  expect_equal(res_bracketed, res_toplevel[, sub], tolerance = 1e-13)
  # ... and, under the identity transform, the raw draws themselves.
  # NOTE: tolerance 1e-6 (not tighter) because cmdstanr's
  # constrain_variables(unconstrain_draws) round-trip differs from the
  # fit$draws() readout at the ~1e-7 level (measured: 1.03e-7 max over all
  # 3000 rows here) -- that is cmdstanr round-trip noise, not package error:
  # the package's own identity transform is exact to 4e-16.
  expect_equal(res_bracketed, fit$draws(sub, format = "draws_matrix"),
               tolerance = 1e-6)
})
