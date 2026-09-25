# ---------------------------------------------------------------------------
# Tests for estimator = "ofs" (Shaby, 2014, open-faced sandwich).
#
# The OFS transform applies Omega = (Q^-1 P)^(1/2) DIRECTLY to the centered
# unconstrained draws (transform_ofs -> adjust_posterior_draws with M =
# Omega), without calibrating to the empirical draw covariance C_emp -- in
# contrast to the closed sandwich, whose transform maps C_emp to V exactly.
# Consequences exercised below:
#   - 1D: se_ofs^2 = Omega^2 * C_emp = (P/Q) * C_emp, while the sandwich
#     gives se_sand^2 = V = P/Q^2; their ratio is C_emp * Q, which is ~ 1
#     for a TRUE likelihood (truelik) where C_emp ~= Q^-1.
#   - bread = "mcmc": Q = solve(C_emp), so in 1D OFS and the closed
#     sandwich coincide exactly (both C_emp^2 * P).
#   - multivariate: the adjusted UNCONSTRAINED draws have covariance
#     Omega C_emp Omega'.
#
# All tests are Stan-backed (gated by skip_stan(); see helper-fit.R).
# truelik has a single real parameter ("theta"), so constrained ==
# unconstrained, and the table SE (sd of the adjusted draws) IS the
# adjusted-draw variance.
# ---------------------------------------------------------------------------

test_that("estimator='ofs' expands the bread x meat x ci grid with ofs.* labels", {
  skip_stan()
  fit <- get_truelik_fit()
  pp <- adjust_pseudo_posterior(fit, "theta",
                            estimator = "ofs",
                            bread = c("mcmc", "hessian_mle"),
                            meat = "score",
                            ci = "quantile",
                            model = get_truelik_model(),
                            data = pp_truelik_dat(),
                            init = list(list(theta = 0.0)),
                            verbose = FALSE)
  expect_setequal(pp$table$method,
                  c("ofs.mcmc.score.quant", "ofs.mle.score.quant"))
})

test_that("ofs runs cleanly with hessian_mle bread and score meat", {
  skip_stan()
  fit <- get_truelik_fit()
  pp <- adjust_pseudo_posterior(fit, "theta",
                                estimator = "ofs",
                                bread = "hessian_mle", meat = "score",
                                ci = "quantile",
                                model = get_truelik_model(),
                                data = pp_truelik_dat(),
                                init = list(list(theta = 0.0)),
                                verbose = FALSE)
  expect_identical(pp$meta$warnings, character(0))
  r <- pp$table[pp$table$method == "ofs.mle.score.quant", ]
  expect_true(is.finite(r$se))
  expect_true(r$se > 0)
})

test_that("ofs 1D: se_ofs^2 / se_sand^2 = C_emp * Q (deviates from 1 by the one-bread gap)", {
  skip_stan()
  fit <- get_truelik_fit()
  mod <- get_truelik_model()
  dat <- pp_truelik_dat()
  init <- list(list(theta = 0.0))

  # OFS
  pp_ofs <- adjust_pseudo_posterior(fit, "theta",
                                    estimator = "ofs", bread = "hessian_mle",
                                    meat = "score", ci = "quantile",
                                    model = mod, data = dat, init = init,
                                    verbose = FALSE)
  se_ofs <- pp_ofs$table$se[pp_ofs$table$method == "ofs.mle.score.quant"]

  # Sandwich (same bread/meat)
  pp_sand <- adjust_pseudo_posterior(fit, "theta",
                                     estimator = "sandwich", bread = "hessian_mle",
                                     meat = "score", ci = "quantile",
                                     model = mod, data = dat, init = init,
                                     verbose = FALSE)
  se_sand <- pp_sand$table$se[pp_sand$table$method == "sandwich.mle.score.quant"]

  # C_emp * Q: compute directly. (The adjust_pseudo_posterior calls above
  # already initialized the fit's model methods via ensure_model_methods.)
  u <- fit$unconstrain_draws(format = "draws_matrix")
  C_emp <- as.numeric(cov(u))
  opt <- mod$optimize(data = utils::modifyList(dat, list(magnitude_adj = 1.0,
                                                         use_priors = 0L)),
                      init = init, jacobian = FALSE, refresh = 0,
                      show_messages = FALSE)
  tryCatch(opt$init_model_methods(), error = function(e) NULL)
  Q <- as.numeric(pseudopost:::get_h_unc(opt, pseudopost:::get_unc(opt),
                                         jacobian = FALSE))

  expected_ratio <- C_emp * Q
  observed_ratio <- se_ofs^2 / se_sand^2
  expect_equal(observed_ratio, expected_ratio, tolerance = 1e-4)
  # For a true likelihood (truelik), C_emp ~ Q^{-1}, so ratio ~ 1.
  expect_true(abs(observed_ratio - 1) < 0.3)
})

test_that("ofs multivariate: adjusted unconstrained draws have cov Omega C_emp Omega'", {
  skip_stan()
  fit <- get_tinypair_fit()
  mod <- get_tinypair_model()
  dat <- pp_make_tinypair_stan_dat()
  init <- pp_tinypair_init()

  # Get Q and P from the MLE point
  opt <- mod$optimize(data = c(dat, list(magnitude_adj = 1.0, use_priors = 0L)),
                      init = init, jacobian = FALSE, refresh = 0,
                      show_messages = FALSE)
  tryCatch(opt$init_model_methods(), error = function(e) NULL)
  tryCatch(fit$init_model_methods(), error = function(e) NULL)

  g <- pseudopost:::pl_ofs(opt, htype = "hessian", jacobian = FALSE)
  Omega <- g$Omega
  Q <- g$Q
  P <- g$P

  # Verify Omega^2 = Q^{-1} P (i.e., Omega %*% Omega = solve(Q) %*% P;
  # the spectral construction satisfies this up to round-off -- the 9x9
  # bread's condition number is ~1e8 here, so allow the same 1e-8 relative
  # slack as the other 9-parameter matrix identities in test-correctness.R).
  expect_equal(Omega %*% Omega, solve(Q) %*% P, tolerance = 1e-8)

  # Apply the OFS transform to the unconstrained draws (no C_emp
  # calibration): the affine map makes cov(u_adj) == Omega C_emp Omega'
  # exactly.
  u <- fit$unconstrain_draws(format = "draws_matrix")
  C_emp <- cov(u)
  u_adj <- pseudopost:::adjust_posterior_draws(u, Omega)
  expect_equal(cov(u_adj), Omega %*% C_emp %*% t(Omega), tolerance = 1e-8,
               ignore_attr = TRUE)
})

test_that("ofs bread='mcmc' in 1D agrees with sandwich (Q = solve(C_emp))", {
  skip_stan()
  fit <- get_truelik_fit()
  mod <- get_truelik_model()
  dat <- pp_truelik_dat()
  init <- list(list(theta = 0.0))

  pp_ofs <- adjust_pseudo_posterior(fit, "theta",
                                    estimator = "ofs", bread = "mcmc",
                                    meat = "score", ci = "quantile",
                                    model = mod, data = dat, init = init,
                                    verbose = FALSE)
  pp_sand <- adjust_pseudo_posterior(fit, "theta",
                                     estimator = "sandwich", bread = "mcmc",
                                     meat = "score", ci = "quantile",
                                     model = mod, data = dat, init = init,
                                     verbose = FALSE)
  expect_identical(pp_ofs$meta$warnings, character(0))
  expect_identical(pp_sand$meta$warnings, character(0))
  se_ofs <- unname(pp_ofs$table$se[pp_ofs$table$method == "ofs.mcmc.score.quant"])
  se_sand <- unname(pp_sand$table$se[pp_sand$table$method == "sandwich.mcmc.score.quant"])
  # Both sides use the identical bread Q = solve(C_emp) and meat P, and in
  # 1D Omega^2 * C_emp == C_emp^2 * P == V, so the SDs of the adjusted draws
  # agree up to round-off.
  expect_equal(se_ofs, se_sand, tolerance = 1e-10)
})

test_that("ofs CI variants reuse the same adjusted draws", {
  skip_stan()
  fit <- get_truelik_fit()
  pp <- adjust_pseudo_posterior(fit, "theta",
                                estimator = "ofs", bread = "mcmc",
                                meat = "score", ci = c("quantile", "wald"),
                                verbose = FALSE)
  expect_identical(pp$meta$warnings, character(0))
  quant <- "ofs.mcmc.score.quant"
  wald <- "ofs.mcmc.score.wald"
  expect_setequal(pp$table$method, c(quant, wald))
  expect_identical(pp$draws[[quant]], pp$draws[[wald]])
})

test_that("ofs bread='hessian_mle' without model/data records a warning", {
  skip_stan()
  fit <- get_truelik_fit()
  pp <- adjust_pseudo_posterior(fit, "theta",
                                estimator = "ofs", bread = "hessian_mle",
                                meat = "score", ci = "quantile",
                                verbose = FALSE)
  expect_false("ofs.mle.score.quant" %in% pp$table$method)
  expect_true(any(grepl("needs a MLE point", pp$meta$warnings)))
})

test_that("ofs validates its axes", {
  skip_stan()
  fit <- get_truelik_fit()
  expect_error(adjust_pseudo_posterior(fit, "theta", estimator = "bogus",
                                       verbose = FALSE), "estimator")
  expect_error(adjust_pseudo_posterior(fit, "theta", estimator = "ofs",
                                       bread = "bogus", verbose = FALSE), "bread")
  expect_error(adjust_pseudo_posterior(fit, "theta", estimator = "ofs",
                                       meat = "bogus", verbose = FALSE), "meat")
})
