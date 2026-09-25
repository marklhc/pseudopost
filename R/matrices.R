# ---------------------------------------------------------------------------
# matrices.R
#
# Low-level linear-algebra and per-fit quantity helpers. These are internal
# building blocks shared by the sandwich (CurvAdj) and IJ (IJSE)
# estimators. The formulas are a direct refactor of the reference
# `pairwise_func.R`; they are preserved exactly.
# ---------------------------------------------------------------------------

#' Matrix (symmetric) square root
#'
#' Computes the symmetric square root of a positive semi-definite matrix via a
#' spectral decomposition, clipping any negligible negative eigenvalues to
#' zero. Internal helper.
#'
#' @param Sigma a numeric matrix.
#' @return a matrix `S` with `S %*% S` equal to the symmetrized `Sigma`.
#' @keywords internal
mat_sqrt <- function(Sigma) {
  Sigma <- (Sigma + t(Sigma)) / 2
  e <- eigen(Sigma, symmetric = TRUE)
  # Ensure symmetric/positive definite by dropping negligible imaginary parts.
  # Note: explicit nrow/ncol because bare diag(<scalar>) is interpreted as an
  # identity-matrix size, which breaks for single-parameter (d = 1) models.
  w <- sqrt(pmax(e$values, 0))
  e$vectors %*% diag(w, nrow = length(w), ncol = length(w)) %*% t(e$vectors)
}

#' Inverse matrix (symmetric) square root
#'
#' Spectral inverse square root, clipping negligible negative eigenvalues to
#' zero (which produces `Inf` on the exact-zero direction, matching the
#' reference behavior). Internal helper.
#'
#' @param Sigma a numeric matrix.
#' @return a matrix `S` with `S %*% Sigma %*% S` equal to the identity on the
#'   positive-eigenvalue subspace.
#' @keywords internal
mat_inv_sqrt <- function(Sigma) {
  Sigma <- (Sigma + t(Sigma)) / 2
  e <- eigen(Sigma, symmetric = TRUE)
  # Explicit nrow/ncol: bare diag(<scalar>) is an identity size, not a
  # length-1 diagonal, which would break for d = 1 models.
  w <- 1 / sqrt(pmax(e$values, 0))
  e$vectors %*% diag(w, nrow = length(w), ncol = length(w)) %*% t(e$vectors)
}

#' Center (posterior-mean / MAP) of a fit in the unconstrained space
#'
#' For a `CmdStanMLE` fit this returns the unconstrained MAP/penalized-MLE; for
#' a `CmdStanMCMC` fit it returns the posterior mean of the unconstrained
#' draws. Internal helper.
#'
#' @param fit a `CmdStanMCMC` or `CmdStanMLE` object.
#' @return a numeric vector of unconstrained parameters.
#' @keywords internal
get_unc <- function(fit) {
  # return MAP/penalized MLE for optimized object, and posterior mean for MCMC
  if (inherits(fit, "CmdStanMLE")) {
    params <- names(fit$variable_skeleton(
      transformed_parameters = FALSE,
      generated_quantities = FALSE
    ))
    map <- lapply(params, fit$mle)
    names(map) <- params
    fit$unconstrain_variables(map) # constrained -> unconstrained vector
  } else {
    colMeans(fit$unconstrain_draws(format = "draws_matrix"))
  }
}

#' Centered per-observation score outer product (the "meat"), unconstrained space
#'
#' Evaluates the Jacobian of the per-observation log-likelihood contributions
#' (the generated quantity `name_lli`) with respect to the unconstrained
#' parameters at `unc`, centers the columns, and returns the outer product.
#' This is the Godambe "meat" `J`. Internal helper.
#'
#' @param fit a `CmdStanMCMC` or `CmdStanMLE` object.
#' @param unc a numeric vector of unconstrained parameters (the center).
#' @param name_lli name of the generated-quantity holding per-observation
#'   log-likelihood contributions.
#' @return a `d x d` matrix (`d` = number of unconstrained parameters).
#' @keywords internal
get_j_unc <- function(fit, unc, name_lli = "log_lik") {
  # plli_unc maps unconstrained params back to constrained via
  # constrain_variables(), then evaluates per-observation log-likelihood
  # contributions.
  plli_unc <- function(u) {
    v <- fit$constrain_variables(u) # unconstrained -> constrained vector
    v[[name_lli]]
  }
  jac <- numDeriv::jacobian(plli_unc, unc, method = "simple")
  # Center the Jacobian by subtracting column means
  jac_c <- sweep(jac, 2, colMeans(jac), "-")
  crossprod(jac_c)
}

#' Posterior-averaged centered per-observation score outer product (the "meat")
#'
#' A draw-based (posterior-average) version of [get_j_unc()]: subsamples
#' `score_draws` posterior draws, evaluates the Jacobian of the per-observation
#' log-likelihood at each, centers the columns, and averages the outer
#' products. Shaby's (2014) draw-averaged moment estimator for the meat.
#' Internal helper.
#'
#' @param fit a `CmdStanMCMC` or `CmdStanMLE` object.
#' @param score_draws number of draws subsampled for the average.
#' @param name_lli generated-quantity name for per-observation
#'   log-likelihood contributions.
#' @return a `d x d` matrix (`d` = number of unconstrained parameters).
#' @keywords internal
get_j_avg <- function(fit, score_draws = 200L, name_lli = "log_lik") {
  u <- fit$unconstrain_draws(format = "draws_matrix")
  s <- min(score_draws, nrow(u))
  idx <- sort(sample.int(nrow(u), size = s, replace = FALSE))
  u_s <- u[idx, , drop = FALSE]
  plli_unc <- function(uvec) fit$constrain_variables(uvec)[[name_lli]]
  J <- matrix(0, nrow = ncol(u), ncol = ncol(u))
  for (i in seq_len(s)) {
    J_t <- numDeriv::jacobian(plli_unc, u_s[i, ], method = "simple")
    J_t <- sweep(J_t, 2, colMeans(J_t), "-")
    J <- J + crossprod(J_t)
  }
  J / s
}

#' Negative Hessian (the "bread"), unconstrained space
#'
#' @param fit a `CmdStanMCMC` or `CmdStanMLE` object.
#' @param unc a numeric vector of unconstrained parameters (the center).
#' @param jacobian logical; pass the Jacobian-adjusted log density to the
#'   Hessian (TRUE) or the raw one (FALSE).
#' @return a `d x d` matrix, the negative Hessian `H`.
#' @keywords internal
get_h_unc <- function(fit, unc, jacobian = TRUE) {
  -fit$hessian(unc, jacobian = jacobian)$hessian
}

#' Empirical posterior variance of the unconstrained draws
#'
#' @param fit a `CmdStanMCMC` object.
#' @return a `d x d` covariance matrix of the unconstrained draws.
#' @keywords internal
get_hinv_postvar <- function(fit) {
  cov(fit$unconstrain_draws(format = "draws_matrix"))
}

#' Jointly compute the meat `J`, bread `H`, and (optionally) `H^-1`
#'
#' Orchestrates the per-quantity helpers. Mirrors the reference `get_jh`,
#' extended with `name_lli` and `jacobian` so the Hessian can be evaluated
#' with or without the Jacobian adjustment. Internal helper.
#'
#' @param fit a `CmdStanMCMC` or `CmdStanMLE` object.
#' @param j logical; compute the meat `J`.
#' @param h logical; compute the bread `H`.
#' @param hinv either logical (`TRUE` to return `H^-1`, `FALSE` to omit) or a
#'   supplied covariance matrix to be inverted into `H`.
#' @param name_lli generated-quantity name for per-observation log-likelihood.
#' @param jacobian logical passed to the Hessian.
#' @return a list with elements `j`, `h`, `hinv` (each a matrix or the
#'   original logical sentinel when not requested).
#' @keywords internal
get_jh <- function(fit, j = TRUE, h = TRUE, hinv = FALSE,
                   name_lli = "log_lik", jacobian = TRUE) {
  if (isTRUE(j) || (isTRUE(h) && is.logical(hinv))) {
    unc <- get_unc(fit)
    if (isTRUE(j)) {
      j <- get_j_unc(fit, unc, name_lli = name_lli)
    }
    if (isTRUE(h) && is.logical(hinv)) {
      h <- get_h_unc(fit, unc, jacobian = jacobian)
    }
  }
  if (!is.logical(hinv) && isTRUE(h)) {
    h <- solve(hinv)
  }
  if (isTRUE(hinv)) {
    hinv <- solve(h)
  }
  list(j = j, h = h, hinv = hinv)
}
