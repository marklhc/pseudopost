# ---------------------------------------------------------------------------
# Pure-R unit tests for the matrix / affine-transform core.
# NO Stan, NO gating -- always run. Complements (does not duplicate) the
# checks in test-basics.R with general-matrix and edge-case coverage.
# ---------------------------------------------------------------------------

test_that("mat_sqrt is a genuine symmetric square root (d = 4)", {
  set.seed(11)
  a <- matrix(rnorm(16), 4)
  s <- crossprod(a) + 4 * diag(4)
  s_sqrt <- pseudopost:::mat_sqrt(s)
  # The square of the sqrt reproduces S.
  expect_equal(s_sqrt %*% s_sqrt, s, tolerance = 1e-10)
  # The result is symmetric.
  expect_equal(s_sqrt, t(s_sqrt), tolerance = 1e-10)
  # Nested sqrt is a fourth root: mat_sqrt(mat_sqrt(S))^4 = S.
  s_fourth <- pseudopost:::mat_sqrt(pseudopost:::mat_sqrt(s))
  expect_equal(s_fourth %*% s_fourth %*% s_fourth %*% s_fourth, s,
               tolerance = 1e-8)
})

test_that("mat_sqrt handles the identity and the d = 1 edge case", {
  expect_equal(pseudopost:::mat_sqrt(diag(1, 3)), diag(3), tolerance = 1e-12)
  # d = 1: a 1x1 matrix holding the scalar sqrt (a bare diag(<scalar>) would
  # be misread as an identity size -- the explicit nrow/ncol in matrices.R is
  # what makes this work).
  expect_equal(pseudopost:::mat_sqrt(matrix(2.5)), matrix(sqrt(2.5)),
               tolerance = 1e-12)
})

test_that("mat_inv_sqrt is the inverse of mat_sqrt (d = 4 and d = 1)", {
  set.seed(12)
  a <- matrix(rnorm(16), 4)
  s <- crossprod(a) + 4 * diag(4)
  expect_equal(pseudopost:::mat_sqrt(s) %*% pseudopost:::mat_inv_sqrt(s),
               diag(4), tolerance = 1e-10)
  expect_equal(pseudopost:::mat_inv_sqrt(matrix(4)), matrix(1 / 2),
               tolerance = 1e-12)
})

test_that("mat_sqrt clips tiny negative eigenvalues instead of erroring", {
  s_psd <- diag(c(1, -1e-12))
  expect_no_error(s_psd_sqrt <- pseudopost:::mat_sqrt(s_psd))
  # The negative-eigenvalue direction is clipped to 0.
  expect_equal(s_psd_sqrt, diag(c(1, 0)), tolerance = 1e-10)
})

test_that("adjust_posterior_draws maps covariance C -> M C M' for arbitrary M", {
  set.seed(13)
  draws_mat <- matrix(rnorm(1000 * 3), 1000, 3)
  colnames(draws_mat) <- c("a", "b", "c")
  # A general (non-symmetric, non-diagonal) transformation.
  m <- matrix(c(1.3, -0.4, 0.2, 0.1, 0.8, 0.5, -0.3, 0.6, 1.1), 3, 3)
  d_adj <- pseudopost::adjust_posterior_draws(draws_mat, m)
  expect_equal(colnames(d_adj), colnames(draws_mat))
  # The affine map acts on centered draws, so the sample covariance is
  # transformed exactly (up to floating point).
  strip_names <- function(x) { dimnames(x) <- NULL; x }
  expect_equal(strip_names(cov(d_adj)),
               strip_names(m %*% cov(draws_mat) %*% t(m)),
               tolerance = 1e-10)
  # The center is preserved.
  expect_equal(colMeans(d_adj), colMeans(draws_mat), tolerance = 1e-10)
})

test_that("ofs_square_root: Omega Omega = Q^-1 P and gives the Godambe target under one-bread", {
  # Two-parameter NON-COMMUTING example (P != Q, and Q and P do not commute):
  # this is the regime where the square-root ORDERING matters. A wrong order
  # (Q^(1/2) S^(1/2) Q^-1/2) squares to P Q^-1 and fails both assertions.
  Q <- matrix(c(3.0, 0.7, 0.7, 2.0), 2)
  P <- matrix(c(1.5, 0.9, 0.9, 2.5), 2)
  stopifnot(max(abs(Q %*% P - P %*% Q)) > 1e-8)
  Omega <- pseudopost:::ofs_square_root(Q, P)
  # Square-root identity: Omega is the (left) square root of Q^-1 P.
  expect_equal(Omega %*% Omega, solve(Q) %*% P, tolerance = 1e-10)
  # One-bread: when the draws' covariance is C_emp = Q^-1, the adjusted-draw
  # covariance Omega C_emp Omega' equals the full Godambe target V =
  # Q^-1 P Q^-1 EXACTLY (the multivariate case the ordering fix targets).
  C_emp <- solve(Q)
  V <- solve(Q) %*% P %*% solve(Q)
  expect_equal(Omega %*% C_emp %*% t(Omega), V, tolerance = 1e-10)
  # Contrast: the transposed (wrong-order) square root does NOT hit V.
  S <- pseudopost:::mat_inv_sqrt(Q) %*% P %*% pseudopost:::mat_inv_sqrt(Q)
  Omega_wrong <- pseudopost:::mat_sqrt(Q) %*% pseudopost:::mat_sqrt(S) %*%
    pseudopost:::mat_inv_sqrt(Q)
  expect_true(max(abs(Omega_wrong %*% C_emp %*% t(Omega_wrong) - V)) > 1e-3)
})

test_that("draw-covariance transform maps the empirical variance to the target (review example)", {
  # Review example: draw var 4, target var 2, inverse-Hessian var 1.
  set.seed(99)
  n <- 10000
  u <- matrix(rnorm(n, 0, 2), n, 1)
  colnames(u) <- "x"
  C_emp <- cov(u)
  T <- matrix(2)
  Hinv <- matrix(1)
  # The corrected map is built from the empirical draw covariance, so the
  # adjusted draws have the TARGET variance (2).
  M_new <- pseudopost:::mat_sqrt(T) %*% pseudopost:::mat_inv_sqrt(C_emp)
  u_new <- pseudopost:::adjust_posterior_draws(u, M_new)
  expect_equal(as.numeric(var(u_new)), 2, tolerance = 1e-3)
  # The OLD map was built from the point's inverse-Hessian (var 1 here); with
  # draw var 4 it produced var 8, not the target 2. Documents the old bug.
  M_old <- pseudopost:::mat_sqrt(T) %*% pseudopost:::mat_inv_sqrt(Hinv)
  u_old <- pseudopost:::adjust_posterior_draws(u, M_old)
  expect_true(abs(as.numeric(var(u_old)) - 8) < 0.2)
})
