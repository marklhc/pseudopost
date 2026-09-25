#' Analytic point Godambe sandwich covariance
#'
#' Computes the Godambe (sandwich) covariance `V = H^-1 J H^-1` in the
#' **unconstrained** parameter space. The "bread" `H` is either the negative
#' Hessian at the evaluation point (`htype = "hessian"`) or the inverse of the
#' empirical posterior draw covariance (`htype = "postvar"`); the "meat" `J`
#' is the centered per-observation score outer product at the evaluation point
#' unless supplied via `j`.
#'
#' @param fit a `CmdStanMCMC` or `CmdStanMLE` object (the score evaluation
#'   point is the unconstrained posterior mean for MCMC, the optimum for MLE).
#' @param htype `"hessian"` (default; use the Hessian at the point) or
#'   `"postvar"` (use the empirical draw covariance as `H^-1`).
#' @param name_lli generated-quantity name for per-observation
#'   log-likelihood contributions.
#' @param jacobian logical; evaluate the Hessian with (TRUE) or without
#'   (FALSE) the Jacobian adjustment. Ignored when `htype = "postvar"`.
#' @param j optional pre-computed meat `J` (skips its computation when given).
#' @param h optional pre-computed bread `H` (skips its computation when given).
#' @return a list with:
#'   \itemize{
#'     \item `V` -- the Godambe covariance `H^-1 J H^-1` (`d x d`).
#'     \item `H` -- the bread (negative Hessian, or `solve` of the postvar).
#'     \item `J` -- the meat.
#'     \item `c_naive` -- the naive covariance at the point: the empirical
#'       posterior variance for `htype = "postvar"`, `H^-1` otherwise.
#'   }
#' @references
#' Godambe, P. C. (1960). An Exact Likelihood Ratio Test for Family of
#'   Distributions. The Annals of Mathematical Statistics, 31(2), 343-349.
#' Godambe, P. C. (1985). Optimal Robust Inference. Biometrika, 72(3),
#'   611-620.
#' @examplesIf identical(Sys.getenv("PP_EXAMPLES_STAN"), "1") && requireNamespace("cmdstanr", quietly = TRUE)
#' library(cmdstanr)
#' model <- cmdstan_model(
#'   system.file("test-models", "truelik.stan", package = "pseudopost"),
#'   force_recompile = TRUE, quiet = TRUE)
#' set.seed(42)
#' dat <- list(N = 100L, y = rnorm(100, 1, 1),
#'             magnitude_adj = 1.0, use_priors = 1L)
#' fit <- model$sample(data = dat, chains = 2L, iter_warmup = 200L,
#'                     iter_sampling = 200L, seed = 42L, refresh = 0,
#'                     show_messages = FALSE)
#' g <- pl_godambe(fit, htype = "hessian")
#' dim(g$V); dim(g$J)
#' @export
pl_godambe <- function(fit, htype = c("hessian", "postvar"),
                       name_lli = "log_lik", jacobian = TRUE,
                       j = NULL, h = NULL) {
  htype <- match.arg(htype)

  if (!is.null(h)) {
    hinv <- solve(h)
  } else if (htype == "postvar") {
    hinv <- get_hinv_postvar(fit)
    h <- solve(hinv)
  } else {
    h <- get_h_unc(fit, get_unc(fit), jacobian = jacobian)
    hinv <- solve(h)
  }

  if (is.null(j)) {
    j <- get_j_unc(fit, get_unc(fit), name_lli = name_lli)
  }

  V <- hinv %*% j %*% hinv
  list(V = V, H = h, J = j, c_naive = hinv)
}

#' Open-Faced Sandwich square root (Shaby, 2014)
#'
#' Computes the open-faced sandwich square root `Omega = (Q^-1 P)^(1/2)` in
#' the **unconstrained** parameter space, where the "bread" `Q` is either the
#' negative Hessian at the evaluation point (`htype = "hessian"`) or the
#' inverse of the empirical posterior draw covariance (`htype = "postvar"`)
#' and the "meat" `P` is the centered per-observation score outer product at
#' the evaluation point unless supplied via `j`.
#'
#' @param fit a `CmdStanMCMC` or `CmdStanMLE` object (the score evaluation
#'   point is the unconstrained posterior mean for MCMC, the optimum for MLE).
#' @param htype `"hessian"` (default; use the Hessian at the point) or
#'   `"postvar"` (use the empirical draw covariance as `Q^-1`).
#' @param name_lli generated-quantity name for per-observation
#'   log-likelihood contributions.
#' @param jacobian logical; evaluate the Hessian with (TRUE) or without
#'   (FALSE) the Jacobian adjustment. Ignored when `htype = "postvar"`.
#' @param j optional pre-computed meat `P` (skips its computation when given).
#' @param h optional pre-computed bread `Q` (skips its computation when given).
#' @return a list with:
#'   \itemize{
#'     \item `Omega` -- the open-faced sandwich square root `(Q^-1 P)^(1/2)`
#'       (`d x d`); apply it to the centered unconstrained draws (as
#'       [adjust_posterior_draws()] does with `M = Omega`) for the OFS
#'       adjustment.
#'     \item `Q` -- the bread (negative Hessian, or `solve` of the postvar).
#'     \item `P` -- the meat.
#'   }
#' @details
#' `Q^-1 P` is not symmetric, but it is similar to the symmetric positive
#' semi-definite matrix `S = Q^-1/2 P Q^-1/2`, so the square root is built
#' from spectral square roots: `Omega = Q^-1/2 S^(1/2) Q^(1/2)` (this
#' satisfies `Omega %*% Omega = Q^-1 P`; the ordering matters -- the similar
#' square root `Q^(1/2) S^(1/2) Q^-1/2` squares to `P Q^-1` instead and does
#' NOT reproduce the Godambe target below).
#'
#' In contrast to [pl_godambe()], which returns the closed Godambe covariance
#' `V = Q^-1 P Q^-1` (and [transform_to_target()] calibrates the draws to that
#' full target exactly), `pl_ofs()` returns the open-faced "one-bread" square
#' root of Shaby (2014): `Omega` is applied directly to the centered draws
#' without any calibration to the empirical draw covariance `C_emp`. The
#' adjusted draws then have covariance `Omega C_emp Omega'`, which under the
#' "one-bread" assumption `C_emp = Q^-1` equals `V` exactly (for any number of
#' parameters); it departs from `V` only to the extent that `C_emp` departs
#' from `Q^-1`.
#' @references
#' Shaby, B. A. (2014). The Open-Faced Sandwich Adjustment for MCMC Using
#'   Estimating Functions. Journal of Computational and Graphical Statistics,
#'   23(3), 853-876. doi:10.1080/10618600.2013.842174.
#' @examplesIf identical(Sys.getenv("PP_EXAMPLES_STAN"), "1") && requireNamespace("cmdstanr", quietly = TRUE)
#' library(cmdstanr)
#' model <- cmdstan_model(
#'   system.file("test-models", "truelik.stan", package = "pseudopost"),
#'   force_recompile = TRUE, quiet = TRUE)
#' set.seed(42)
#' dat <- list(N = 100L, y = rnorm(100, 1, 1),
#'             magnitude_adj = 1.0, use_priors = 1L)
#' fit <- model$sample(data = dat, chains = 2L, iter_warmup = 200L,
#'                     iter_sampling = 200L, seed = 42L, refresh = 0,
#'                     show_messages = FALSE)
#' g <- pl_ofs(fit, htype = "hessian")
#' dim(g$Omega); dim(g$Q); dim(g$P)
#' @export
pl_ofs <- function(fit, htype = c("hessian", "postvar"),
                   name_lli = "log_lik", jacobian = TRUE,
                   j = NULL, h = NULL) {
  htype <- match.arg(htype)

  if (!is.null(h)) {
    Q <- h
  } else if (htype == "postvar") {
    Q <- solve(get_hinv_postvar(fit))
  } else {
    Q <- get_h_unc(fit, get_unc(fit), jacobian = jacobian)
  }

  if (is.null(j)) {
    P <- get_j_unc(fit, get_unc(fit), name_lli = name_lli)
  } else {
    P <- j
  }

  Omega <- ofs_square_root(Q, P)
  list(Omega = Omega, Q = Q, P = P)
}

#' Open-faced sandwich square root from explicit bread and meat
#'
#' Internal building block for [pl_ofs()]: builds `Omega = (Q^-1 P)^(1/2)`
#' from a symmetric positive-definite bread `Q` and a symmetric positive
#' semi-definite meat `P` via spectral square roots,
#' `Omega = Q^-1/2 S^(1/2) Q^(1/2)` with `S = Q^-1/2 P Q^-1/2`. The ordering is
#' the one that satisfies `Omega %*% Omega = Q^-1 P` (so that, under the
#' one-bread assumption `C_emp = Q^-1`, the adjusted-draw covariance
#' `Omega C_emp Omega'` equals the Godambe target `Q^-1 P Q^-1`).
#'
#' @param Q a `d x d` symmetric positive-definite bread matrix.
#' @param P a `d x d` symmetric positive semi-definite meat matrix.
#' @return the `d x d` open-faced sandwich square root `Omega`.
#' @keywords internal
ofs_square_root <- function(Q, P) {
  Q_sqrt <- mat_sqrt(Q)
  Q_inv_sqrt <- mat_inv_sqrt(Q)
  S <- Q_inv_sqrt %*% P %*% Q_inv_sqrt
  Q_inv_sqrt %*% mat_sqrt(S) %*% Q_sqrt
}
