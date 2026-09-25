# ---------------------------------------------------------------------------
# Method identity: the 'proxy' IJ target IS the infinitesimal
# jackknife (IJ). We assert it equals the independent-data formula of the CRAN
# package IJSE (Ji, Lee & Rabe-Hesketh 2024; Giordano & Broderick 2023),
# translated verbatim:
#   inf <- N * cov(draws, post_log)
#   V   <- cov(t(inf)) / N
# and confirm our full pipeline's robust SEs track that formula in the
# constrained parameter space.
#
# Gated by skip_stan(); reuses the session-cached tinypair fit (helper-fit.R).
# ---------------------------------------------------------------------------

test_that("'proxy' IJ target == infinitesimal jackknife (CRAN IJSE formula)", {
  skip_stan()
  fit <- get_tinypair_fit()
  fit$init_model_methods()
  lli <- fit$draws("log_lik", format = "draws_matrix")
  u   <- fit$unconstrain_draws(format = "draws_matrix")
  N   <- ncol(lli)

  # Verbatim IJSE::IJ_se (independent case) robust covariance, applied to the
  # unconstrained draws (same space our 'proxy' target is returned in):
  V_ijse <- cov(t(N * cov(u, lli))) / N

  expect_equal(unname(pl_ij(fit, method = "proxy")),
               unname(V_ijse), tolerance = 1e-10)
})

test_that("IJ robust SEs track the IJSE formula (constrained space)", {
  skip_stan()
  fit <- get_tinypair_fit()
  fit$init_model_methods()
  lli    <- fit$draws("log_lik", format = "draws_matrix")
  N      <- ncol(lli)
  params <- pp_tinypair_param_names()

  # IJSE formula applied directly to the naive constrained draws:
  constr  <- fit$draws(params, format = "draws_matrix")
  se_ijse <- sqrt(pmax(diag(cov(t(N * cov(constr, lli)))) / N, 0))
  names(se_ijse) <- params

  pp      <- adjust_pseudo_posterior(fit, params, methods = "ijse2", verbose = FALSE)
  se_ours <- apply(pp$draws$ijse2, 2, sd)
  names(se_ours) <- params

  # The affine rescale -> constrain pipeline is a first-order approximation
  # of the IJ computed directly in the constrained space, so per-parameter
  # agreement is loose; deviations are driven by the nonlinear delta transforms.
  rel <- abs(se_ours - se_ijse) / se_ijse
  expect_true(all(rel < 0.30),
              info = paste(sprintf("%s: ours=%.4f ijse=%.4f (rel=%.3f)",
                                   names(rel), se_ours, se_ijse, rel),
                           collapse = "; "))
})
