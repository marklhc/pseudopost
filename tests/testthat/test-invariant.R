# ---------------------------------------------------------------------------
# CROWN JEWEL (invariant): when the "pseudo-likelihood" is actually a TRUE
# per-observation likelihood (truelik.stan), the information equality
# B = H holds, so every sandwich / IJ correction must reduce to ~ the
# identity: corrected SE ~= noadj SE, and independent estimators agree.
#
# Gated by skip_stan(); reuses the session-cached truelik fit (helper-fit.R).
# ---------------------------------------------------------------------------

test_that("IJ correction ~ identity for a true likelihood (proxy)", {
  skip_stan()
  fit <- get_truelik_fit()
  noadj_se <- sd(fit$draws("theta", format = "draws_matrix"))
  expect_true(is.finite(noadj_se) && noadj_se > 0)

  pp <- adjust_pseudo_posterior(fit, "theta",
                                methods = c("noadj", "ijse", "ijse2"),
                                ij_method = "proxy", verbose = FALSE)
  t_p <- pp$table
  # Internal consistency: the noadj row reproduces the raw draws summary.
  expect_equal(unname(t_p$se[t_p$method == "noadj"]), noadj_se)
  # truelik has no parameter transforms, so the affine transform (in the
  # unconstrained = constrained space) preserves the mean exactly.
  expect_true(max(abs(t_p$est - t_p$est[1])) < 1e-10,
              info = paste("ests:", paste(t_p$est, collapse = ", ")))
  # ijse and ijse2 SEs must be finite and within 30% of the noadj SE.
  for (m in c("ijse", "ijse2")) {
    se <- unname(t_p$se[t_p$method == m])
    expect_true(is.finite(se), info = paste(m, "se is not finite"))
    expect_true(abs(se - noadj_se) / noadj_se < 0.30,
                info = sprintf("%s se = %.6g not within 30%% of noadj_se = %.6g",
                               m, se, noadj_se))
  }
})

test_that("IJ correction ~ identity for a true likelihood (score)", {
  skip_stan()
  fit <- get_truelik_fit()
  noadj_se <- sd(fit$draws("theta", format = "draws_matrix"))

  # The score variant warns about its finite-difference cost; assert that.
  # NOTE: testthat 3e's expect_warning() returns the captured CONDITION, not
  # the value of the expression, so capture the warning by hand (keeping the
  # result) and muffle it so it does not leak into the test results.
  w <- NULL
  pp <- withCallingHandlers(
    adjust_pseudo_posterior(fit, "theta", methods = c("ijse", "ijse2"),
                            ij_method = "score", score_draws = 60L,
                            verbose = FALSE),
    warning = function(cnd) {
      if (is.null(w)) w <<- conditionMessage(cnd)
      invokeRestart("muffleWarning")
    })
  expect_match(w, "expensive")
  t_s <- pp$table
  for (m in c("ijse", "ijse2")) {
    se <- unname(t_s$se[t_s$method == m])
    expect_true(is.finite(se), info = paste(m, "se is not finite"))
    expect_true(abs(se - noadj_se) / noadj_se < 0.30,
                info = sprintf("%s se = %.6g not within 30%% of noadj_se = %.6g",
                               m, se, noadj_se))
  }

  # The two independent estimators (proxy vs score) must converge on the same
  # corrected SE (within 25%).
  pp_proxy <- adjust_pseudo_posterior(fit, "theta", methods = "ijse2",
                                      ij_method = "proxy", verbose = FALSE)
  se_proxy <- unname(pp_proxy$table$se)
  se_score <- unname(t_s$se[t_s$method == "ijse2"])
  expect_true(abs(se_proxy - se_score) / min(se_proxy, se_score) < 0.25,
              info = sprintf("proxy ijse2 se = %.6g vs score ijse2 se = %.6g (rel. diff > 25%%)",
                             se_proxy, se_score))
})

test_that("curvadj family ~ identity for a true likelihood", {
  skip_stan()
  fit <- get_truelik_fit()
  noadj_se <- sd(fit$draws("theta", format = "draws_matrix"))

  pp <- adjust_pseudo_posterior(fit, "theta",
                                methods = c("curvadj", "curvadj2", "curvadj3"),
                                model = get_truelik_model(),
                                data = pp_truelik_dat(),
                                init = list(list(theta = 0.0)),
                                verbose = FALSE)
  # model/data supplied => no fallback notes.
  expect_length(pp$meta$notes, 0)
  t_c <- pp$table
  expect_setequal(t_c$method, c("curvadj", "curvadj2", "curvadj3"))
  for (m in c("curvadj", "curvadj2", "curvadj3")) {
    se <- unname(t_c$se[t_c$method == m])
    expect_true(is.finite(se), info = paste(m, "se is not finite"))
    expect_true(abs(se - noadj_se) / noadj_se < 0.30,
                info = sprintf("%s se = %.6g not within 30%% of noadj_se = %.6g",
                               m, se, noadj_se))
  }
})

test_that("adjust_pseudo_posterior returns a well-formed pseudo_post (truelik)", {
  skip_stan()
  fit <- get_truelik_fit()
  pp <- adjust_pseudo_posterior(fit, "theta",
                                methods = c("noadj", "ijse2"), verbose = FALSE)
  expect_s3_class(pp, "pseudo_post")
  expect_named(pp, c("table", "draws", "meta"))
  expect_setequal(names(pp$table),
                  c("method", "param", "est", "se", "lo", "hi"))
  expect_named(pp$draws, c("noadj", "ijse2"))
  expect_equal(dim(pp$draws$noadj), c(2000L, 1L))
  expect_equal(pp$meta$ci_level, 0.90)
  expect_equal(pp$meta$ij_method, "proxy")
  expect_identical(pp$meta$warnings, character(0))
  # The 90% CI brackets the estimate.
  expect_true(all(pp$table$lo < pp$table$est))
  expect_true(all(pp$table$est < pp$table$hi))
})
