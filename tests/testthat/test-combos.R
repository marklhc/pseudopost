# ---------------------------------------------------------------------------
# Systematic interface (estimator / bread / meat / ci) of
# adjust_pseudo_posterior(): grid expansion + label scheme, bit-identical
# equivalence with the legacy `methods` presets, bread/meat handling for the
# "ij" estimator, the previously-missing sandwich x wald cell, the score_avg
# meat, input validation, and `methods` precedence.
#
# Gated by skip_stan(); reuses the session-cached truelik fit (helper-fit.R).
# truelik has a single parameter ("theta"), so every table has one row per
# method and SEs are compared on that single row.
# ---------------------------------------------------------------------------

test_that("new mode: a ci vector expands to ij.<ci> labels", {
  skip_stan()
  fit <- get_truelik_fit()
  pp <- adjust_pseudo_posterior(fit, "theta",
                                estimator = "ij",
                                ci = c("quantile", "wald"),
                                verbose = FALSE)
  expect_s3_class(pp, "pseudo_post")
  expect_setequal(pp$table$method, c("ij.quant", "ij.wald"))
  expect_equal(pp$meta$methods, c("ij.quant", "ij.wald"))
})

test_that("new mode: a bread vector expands to the labeled sandwich grid", {
  skip_stan()
  fit <- get_truelik_fit()
  pp <- adjust_pseudo_posterior(fit, "theta",
                                estimator = "sandwich",
                                bread = c("mcmc", "hessian_mle", "hessian_map"),
                                meat = "score",
                                ci = "quantile",
                                model = get_truelik_model(),
                                data = pp_truelik_dat(),
                                init = list(list(theta = 0.0)),
                                verbose = FALSE)
  expect_setequal(pp$table$method,
                  c("sandwich.mcmc.score.quant",
                    "sandwich.mle.score.quant",
                    "sandwich.map.score.quant"))
})

test_that("new interface is bit-identical to the legacy presets (truelik)", {
  skip_stan()
  fit <- get_truelik_fit()
  mod <- get_truelik_model()
  dat <- pp_truelik_dat()
  init <- list(list(theta = 0.0))

  pp_new <- function(estimator, ci, bread = "mcmc") {
    adjust_pseudo_posterior(fit, "theta",
                            estimator = estimator, bread = bread, ci = ci,
                            model = mod, data = dat, init = init,
                            verbose = FALSE)
  }

  cases <- list(
    list(est = "ij", bread = "mcmc", ci = "quantile",
         label = "ij.quant", legacy = "ijse2"),
    list(est = "ij", bread = "mcmc", ci = "wald",
         label = "ij.wald", legacy = "ijse"),
    list(est = "sandwich", bread = "mcmc", ci = "quantile",
         label = "sandwich.mcmc.score.quant", legacy = "curvadj3"),
    list(est = "sandwich", bread = "hessian_mle", ci = "quantile",
         label = "sandwich.mle.score.quant", legacy = "curvadj"),
    list(est = "sandwich", bread = "hessian_map", ci = "quantile",
         label = "sandwich.map.score.quant", legacy = "curvadj2")
  )

  for (cs in cases) {
    pp_n <- pp_new(cs$est, cs$ci, cs$bread)
    pp_l <- adjust_pseudo_posterior(fit, "theta", methods = cs$legacy,
                                    model = mod, data = dat, init = init,
                                    verbose = FALSE)
    # Neither side may have silently failed a method.
    expect_identical(pp_n$meta$warnings, character(0))
    expect_identical(pp_l$meta$warnings, character(0))
    se_n <- unname(pp_n$table$se[pp_n$table$method == cs$label])
    se_l <- unname(pp_l$table$se[pp_l$table$method == cs$legacy])
    expect_length(se_n, 1)
    expect_length(se_l, 1)
    expect_equal(se_n, se_l, tolerance = 1e-12,
                 info = sprintf("new '%s' vs legacy '%s': %.17g vs %.17g",
                                cs$label, cs$legacy, se_n, se_l))
  }
})

test_that("estimator='ij' ignores bread and meat", {
  skip_stan()
  fit <- get_truelik_fit()
  pp_def <- adjust_pseudo_posterior(fit, "theta", estimator = "ij",
                                    verbose = FALSE)
  pp_bm <- adjust_pseudo_posterior(fit, "theta", estimator = "ij",
                                   bread = "hessian_map", meat = "score_avg",
                                   verbose = FALSE)
  # bread/meat are placeholders for "ij": the label does not include them.
  expect_setequal(pp_bm$table$method, "ij.quant")
  se_def <- unname(pp_def$table$se[pp_def$table$method == "ij.quant"])
  se_bm <- unname(pp_bm$table$se[pp_bm$table$method == "ij.quant"])
  expect_equal(se_bm, se_def, tolerance = 1e-12)
})

test_that("the sandwich x wald cell runs and reports a Wald interval", {
  skip_stan()
  fit <- get_truelik_fit()
  pp <- adjust_pseudo_posterior(fit, "theta",
                                estimator = "sandwich", bread = "mcmc",
                                ci = "wald", ci_level = 0.90,
                                model = get_truelik_model(),
                                data = pp_truelik_dat(),
                                init = list(list(theta = 0.0)),
                                verbose = FALSE)
  expect_setequal(pp$table$method, "sandwich.mcmc.score.wald")
  r <- pp$table
  expect_true(is.finite(r$se), info = "sandwich wald se is not finite")
  # Wald CI: hi - lo must be exactly 2 * z * se at ci_level = 0.90.
  expect_true(abs((r$hi - r$lo) - 2 * qnorm(0.95) * r$se) < 1e-8,
              info = sprintf("hi - lo = %.12g vs 2*qnorm(0.95)*se = %.12g",
                             r$hi - r$lo, 2 * qnorm(0.95) * r$se))
})

test_that("sandwich meat='score_avg' runs with the scoreavg label", {
  skip_stan()
  fit <- get_truelik_fit()
  # The new-mode score_avg path calls get_j_avg() directly and emits no
  # warning (unlike legacy ij_method = "score"); assert both.
  pp <- adjust_pseudo_posterior(fit, "theta",
                                estimator = "sandwich", bread = "mcmc",
                                meat = "score_avg", ci = "quantile",
                                score_draws = 30L,
                                model = get_truelik_model(),
                                data = pp_truelik_dat(),
                                init = list(list(theta = 0.0)),
                                verbose = FALSE)
  expect_setequal(pp$table$method, "sandwich.mcmc.scoreavg.quant")
  expect_identical(pp$meta$warnings, character(0))
  expect_true(is.finite(pp$table$se), info = "score_avg se is not finite")
})

test_that("CI variants reuse the same adjusted draws and score sample", {
  skip_stan()
  fit <- get_truelik_fit()
  set.seed(123)
  pp <- adjust_pseudo_posterior(
    fit, "theta", estimator = "sandwich", bread = "mcmc",
    meat = "score_avg", ci = c("quantile", "wald"), score_draws = 8L,
    verbose = FALSE)
  expect_identical(pp$meta$warnings, character(0))
  quant <- "sandwich.mcmc.scoreavg.quant"
  wald <- "sandwich.mcmc.scoreavg.wald"
  expect_setequal(pp$table$method, c(quant, wald))
  expect_identical(pp$draws[[quant]], pp$draws[[wald]])
  expect_equal(pp$table$se[pp$table$method == quant],
               pp$table$se[pp$table$method == wald])
})

test_that("the new interface validates its axes", {
  skip_stan()
  fit <- get_truelik_fit()
  expect_error(adjust_pseudo_posterior(fit, "theta", estimator = "bogus",
                                       verbose = FALSE),
               "estimator")
  expect_error(adjust_pseudo_posterior(fit, "theta", estimator = "sandwich",
                                       bread = "bogus", verbose = FALSE),
               "bread")
  # bread is ALWAYS validated, even when estimator = "ij" ignores it.
  expect_error(adjust_pseudo_posterior(fit, "theta", estimator = "ij",
                                       bread = "bogus", verbose = FALSE),
               "bread")
  expect_error(adjust_pseudo_posterior(fit, "theta", ci = "bogus",
                                        verbose = FALSE),
               "ci")
  expect_error(adjust_pseudo_posterior(fit, "theta", meat = "bogus",
                                        verbose = FALSE),
               "meat")
})

test_that("methods= takes precedence over the new interface", {
  skip_stan()
  fit <- get_truelik_fit()
  pp_leg <- adjust_pseudo_posterior(fit, "theta", methods = "ijse2",
                                    verbose = FALSE)
  # The legacy preset label, NOT the systematic label.
  expect_setequal(pp_leg$table$method, "ijse2")
  expect_setequal(pp_leg$meta$methods, "ijse2")
  pp_new <- adjust_pseudo_posterior(fit, "theta", estimator = "ij",
                                    ci = "quantile", verbose = FALSE)
  se_leg <- unname(pp_leg$table$se[pp_leg$table$method == "ijse2"])
  se_new <- unname(pp_new$table$se[pp_new$table$method == "ij.quant"])
  expect_equal(se_leg, se_new, tolerance = 1e-12)
})

test_that("methods='noadj' reproduces the raw draws summary", {
  skip_stan()
  fit <- get_truelik_fit()
  noadj_se <- sd(fit$draws("theta", format = "draws_matrix"))
  pp <- adjust_pseudo_posterior(fit, "theta", methods = "noadj",
                                verbose = FALSE)
  expect_setequal(pp$table$method, "noadj")
  expect_equal(unname(pp$table$se[pp$table$method == "noadj"]), noadj_se)
})
