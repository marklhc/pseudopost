# ---------------------------------------------------------------------------
# Correctness tests for the corrected hessian-bread (CurvAdj) transform.
#
# The corrected transform_to_target() builds its map from the empirical draw
# covariance C_emp = cov(fit$unconstrain_draws()) -- NOT from the point's
# inverse-Hessian -- so the adjusted UNCONSTRAINED draws have covariance
# exactly equal to the Godambe target V. These tests pin that behavior down
# (test-regression.R deliberately no longer compares curvadj/curvadj2 against
# the old H^-1-based reference rows).
#
# All tests are Stan-backed (gated by skip_stan(); see helper-fit.R).
# ---------------------------------------------------------------------------

# truelik has a single real parameter (theta), so constrained == unconstrained:
# the table SE^2 (sd of the adjusted draws) IS the adjusted-draw variance.

test_that("truelik (MLE): curvadj adjusted-draw variance equals the Godambe target", {
  skip_stan()
  mod <- get_truelik_model()
  dat <- pp_truelik_dat()
  fit <- get_truelik_fit()

  opt_mle <- mod$optimize(
    data = utils::modifyList(dat, list(magnitude_adj = 1.0, use_priors = 0L)),
    init = list(list(theta = 0.0)), jacobian = FALSE,
    refresh = 0, show_messages = FALSE)
  tryCatch(opt_mle$init_model_methods(), error = function(e) NULL)
  V_mle <- as.numeric(pseudopost:::pl_godambe(opt_mle, htype = "hessian",
                                              name_lli = "log_lik",
                                              jacobian = FALSE)$V)

  pp <- adjust_pseudo_posterior(fit, "theta", methods = "curvadj",
                                model = mod, data = dat,
                                init = list(list(theta = 0.0)),
                                verbose = FALSE)
  expect_identical(pp$meta$warnings, character(0))
  se <- pp$table$se[pp$table$method == "curvadj"]
  # KEY assertion: the corrected map is built from the draw covariance, so the
  # adjusted-draw variance (SE^2) equals the Godambe target exactly.
  expect_equal(se^2, V_mle, tolerance = 1e-6)
  # Sanity anchor: truelik is a TRUE likelihood, so the corrected SE must be
  # close to the unadjusted posterior SE.
  se_noadj <- sd(fit$draws("theta", format = "draws_matrix"))
  expect_true(abs(se / se_noadj - 1) < 0.3)
})

test_that("truelik (MAP): curvadj2 adjusted-draw variance equals the Godambe target", {
  skip_stan()
  mod <- get_truelik_model()
  dat <- pp_truelik_dat()
  fit <- get_truelik_fit()

  opt_map <- mod$optimize(
    data = utils::modifyList(dat, list(magnitude_adj = 1.0, use_priors = 1L)),
    init = list(list(theta = 0.0)), jacobian = TRUE,
    refresh = 0, show_messages = FALSE)
  tryCatch(opt_map$init_model_methods(), error = function(e) NULL)
  V_map <- as.numeric(pseudopost:::pl_godambe(opt_map, htype = "hessian",
                                              name_lli = "log_lik",
                                              jacobian = TRUE)$V)

  pp <- adjust_pseudo_posterior(fit, "theta", methods = "curvadj2",
                                model = mod, data = dat,
                                init = list(list(theta = 0.0)),
                                verbose = FALSE)
  expect_identical(pp$meta$warnings, character(0))
  se <- pp$table$se[pp$table$method == "curvadj2"]
  # KEY assertion (as above, at the MAP with the Jacobian adjustment).
  expect_equal(se^2, V_map, tolerance = 1e-6)
  # Sanity anchor: truelik is a TRUE likelihood.
  se_noadj <- sd(fit$draws("theta", format = "draws_matrix"))
  expect_true(abs(se / se_noadj - 1) < 0.3)
})

test_that("tinypair (MLE): corrected unconstrained draw transform has covariance == V", {
  skip_stan()
  fit <- get_tinypair_fit()
  mod <- get_tinypair_model()
  dat <- pp_make_tinypair_stan_dat()
  init <- pp_tinypair_init()
  opt_mle <- mod$optimize(
    data = c(dat, list(magnitude_adj = 1.0, use_priors = 0L)),
    init = init, jacobian = FALSE,
    refresh = 0, show_messages = FALSE)

  tryCatch(fit$init_model_methods(), error = function(e) NULL)
  tryCatch(opt_mle$init_model_methods(), error = function(e) NULL)

  u <- fit$unconstrain_draws(format = "draws_matrix")
  C_emp <- cov(u)
  g <- pseudopost:::pl_godambe(opt_mle, htype = "hessian",
                               name_lli = "log_lik", jacobian = FALSE)
  V <- g$V

  # The corrected map: build from the empirical draw covariance.
  M <- pseudopost:::mat_sqrt(V) %*% pseudopost:::mat_inv_sqrt(C_emp)
  u_adj <- pseudopost:::adjust_posterior_draws(u, M)
  expect_equal(cov(u_adj), V, tolerance = 1e-8, ignore_attr = TRUE)

  # Contrast: the OLD map (built from the point's inverse-Hessian) is far off
  # for this pseudo-likelihood model (C_emp != H^-1 by ~100%).
  M_old <- pseudopost:::mat_sqrt(V) %*% pseudopost:::mat_inv_sqrt(g$c_naive)
  expect_true(max(abs(cov(pseudopost:::adjust_posterior_draws(u, M_old)) - V)) >
                0.5 * max(abs(diag(V))))
})

test_that("curvadj/curvadj2 require an explicit MLE/MAP point (no posterior-mean fallback)", {
  skip_stan()
  fit <- get_truelik_fit()

  # Without model/data/opt: the hessian-bread methods FAIL (recorded in
  # pp$meta$warnings); the MCMC posterior mean is not substituted.
  pp <- adjust_pseudo_posterior(fit, "theta", methods = "curvadj",
                                verbose = FALSE)
  expect_false("curvadj" %in% pp$table$method)
  expect_true(any(grepl("needs a MLE point", pp$meta$warnings)))

  pp2 <- adjust_pseudo_posterior(fit, "theta", methods = "curvadj2",
                                 verbose = FALSE)
  expect_false("curvadj2" %in% pp2$table$method)
  expect_true(any(grepl("needs a MAP point", pp2$meta$warnings)))

  # With model/data it succeeds cleanly.
  pp3 <- adjust_pseudo_posterior(fit, "theta", methods = "curvadj",
                                 model = get_truelik_model(),
                                 data = pp_truelik_dat(),
                                 init = list(list(theta = 0.0)),
                                 verbose = FALSE)
  expect_identical(pp3$meta$warnings, character(0))
  expect_true("curvadj" %in% pp3$table$method)
})
