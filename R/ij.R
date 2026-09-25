#' Infinitesimal jackknife (IJ) target covariance
#'
#' Estimates a target covariance in the **unconstrained** parameter space
#' (`d` = number of unconstrained parameters). `method = "proxy"` computes the
#' draw-based infinitesimal jackknife; `method = "score"` uses a score-averaged
#' meat and a Hessian bread.
#'
#' @param fit a `CmdStanMCMC` object.
#' @param name_lli name of the generated quantity holding per-observation
#'   log-likelihood contributions (default `"log_lik"`).
#' @param method `"proxy"` (default; fast regression proxy) or `"score"`
#'   (literal per-observation scores, expensive).
#' @param score_draws number of draws subsampled for the `"score"` variant
#'   (capped at the number of available draws).
#' @return a `d x d` target covariance matrix in the unconstrained space.
#' @details
#' `method = "proxy"` computes the infinitesimal jackknife (IJ) covariance:
#' form `n_obs * cov(log_lik, unconstrained_draws)` across posterior draws, then
#' take its covariance across observations divided by `n_obs`. This is the
#' draw-based influence-function approach of Giordano & Broderick (2023).
#' `method = "score"` averages per-observation score outer products across
#' draws, then sandwiches that meat with the inverse posterior Hessian.
#' @references
#' Shaby, B. A. (2014). The Open-Faced Sandwich Adjustment for MCMC Using
#'   Estimating Functions. Journal of Computational and Graphical Statistics,
#'   23(3), 853-876. doi:10.1080/10618600.2013.842174.
#' Giordano, R., & Broderick, T. (2023). The Bayesian Infinitesimal Jackknife
#'   for Variance. arXiv preprint arXiv:2305.06466.
#'   doi:10.48550/arXiv.2305.06466.
#' Ji, F., Lee, J., & Rabe-Hesketh, S. (2024). Valid standard errors for
#'   Bayesian quantile regression with clustered and independent data.
#'   arXiv preprint arXiv:2407.09772. doi:10.48550/arXiv.2407.09772.
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
#' B <- pl_ij(fit, method = "proxy")   # fast regression-proxy variant
#' dim(B)
#' @export
pl_ij <- function(fit, name_lli = "log_lik",
                         method = c("proxy", "score"),
                         score_draws = 200L) {
  method <- match.arg(method)

  lli <- fit$draws(name_lli, format = "draws_matrix")
  n_obs <- ncol(lli)
  u <- fit$unconstrain_draws(format = "draws_matrix")

  if (method == "proxy") {
    # EXACT reproduction of the reference cov_ij / cov_ij_direct, applied to
    # the unconstrained draws:
    #   inf_score = n_obs * cov(lli, u)          [m_pseudo x d]
    #   B         = cov(inf_score) / n_obs       [d x d]
    inf_score <- n_obs * cov(lli, u)
    B <- cov(inf_score) / n_obs
  } else {
    # Literal Shaby (2014)-style Godambe sandwich, returned as the *target*
    # covariance so it plugs into the caller exactly like the "proxy" variant.
    #
    #   meat   B = posterior-averaged per-observation score outer product
    #              ([get_j_avg()], Shaby's draw-averaged moment estimator)
    #   bread  A = posterior Hessian (negative) at the center
    #   target T = A^{-1} B A^{-1}
    #
    # For a TRUE likelihood B ~ A (information equality), so T ~ A^{-1} ~ the
    # naive posterior covariance and the affine transform is ~ the identity.
    warning("pl_ij(method='score') is expensive: it performs ~ ",
            min(score_draws, nrow(u)), " x (", ncol(u) + 1L, ") model ",
            "evaluations via finite differences. This may take a long time.",
            call. = FALSE)
    B <- get_j_avg(fit, score_draws = score_draws, name_lli = name_lli)
    A <- get_h_unc(fit, get_unc(fit), jacobian = TRUE)
    Ainv <- solve(A)
    B <- Ainv %*% B %*% Ainv
  }

  dim(B) <- c(ncol(u), ncol(u))
  B
}
