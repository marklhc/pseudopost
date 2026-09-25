# Pure-R core tests; Stan-backed tests live in the other test files.

test_that("mat_sqrt / mat_inv_sqrt are mutual inverses on SPD matrices", {
  set.seed(1)
  a <- matrix(rnorm(16), 4)
  s <- crossprod(a) + diag(4)
  s_sqrt <- pseudopost:::mat_sqrt(s)
  s_inv_sqrt <- pseudopost:::mat_inv_sqrt(s)
  expect_equal(s_sqrt %*% s_sqrt, s, tolerance = 1e-10)
  expect_equal(s_sqrt %*% s_inv_sqrt, diag(4), tolerance = 1e-10)
})

test_that("mat_sqrt / mat_inv_sqrt handle the d = 1 (scalar) edge case", {
  expect_equal(pseudopost:::mat_sqrt(matrix(4)), matrix(2))
  expect_equal(pseudopost:::mat_inv_sqrt(matrix(4)), matrix(1 / 2))
})

test_that("adjust_posterior_draws preserves the center and rescales covariance", {
  set.seed(2)
  d <- matrix(rnorm(2000 * 3), 2000, 3)
  colnames(d) <- c("a", "b", "c")
  c_emp <- cov(d)
  # Target = scalar multiple of the empirical covariance, so the affine map
  # is an exact rescale (the sqrt matrices commute).
  t_target <- 2 * c_emp
  m <- pseudopost:::mat_sqrt(t_target) %*% pseudopost:::mat_inv_sqrt(c_emp)
  d_adj <- pseudopost::adjust_posterior_draws(d, m)
  expect_equal(colnames(d_adj), colnames(d))
  expect_equal(colMeans(d_adj), colMeans(d), tolerance = 1e-10)
  expect_equal(cov(d_adj), t_target, tolerance = 1e-8)
})

test_that("adjust_posterior_draws is the identity for M = I", {
  d <- matrix(rnorm(50), 10, 5)
  expect_equal(pseudopost::adjust_posterior_draws(d, diag(5)), d)
})

make_pp <- function() {
  structure(
    list(
      table = data.frame(
        method = "noadj", param = "a",
        est = 1.23456, se = 0.1, lo = 0.8, hi = 1.2,
        stringsAsFactors = FALSE
      ),
      draws = list(noadj = matrix(1:4, 2, 2)),
      meta = list(fit = NULL, methods = "noadj", ci_level = 0.90,
                  ij_method = "proxy", warnings = character(0),
                  notes = "a note")
    ),
    class = "pseudo_post"
  )
}

test_that("print / summary / as.data.frame work on a pseudo_post object", {
  pp <- make_pp()
  expect_output(print(pp), "Pseudo-posterior uncertainty adjustment")
  expect_output(print(summary(pp)), "pseudo_post summary")
  expect_identical(as.data.frame(pp), pp$table)
  expect_output(print(pp), "a note") # note is surfaced
})

test_that("check_pl_requirements rejects non-cmdstanr objects", {
  expect_error(check_pl_requirements(42), "CmdStanMCMC or CmdStanMLE")
  expect_error(check_pl_requirements(list(a = 1)), "CmdStanMCMC or CmdStanMLE")
})

test_that("pl_ij dispatches on method", {
  expect_error(pl_ij(NULL, method = "bogus"), "should be one of")
})
