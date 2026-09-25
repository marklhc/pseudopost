# ---------------------------------------------------------------------------
# magnitude.R
#
# Scalar magnitude (scale) adjustment. A single scalar `m` rescales the whole
# pseudo-likelihood (`target += magnitude_adj * pll`) to correct the overall
# (global) posterior scale. This is a coarse scalar correction of the overall
# scale; it does not in general reproduce a direction-dependent multivariate
# sandwich covariance. The scalar is the Godambe ratio
#     m = nrow(H) / sum(diag(solve(H, J)))
# evaluated at an (optionally freshly optimized) point. This is a direct
# refactor of the reference `get_mag_adj`.
# ---------------------------------------------------------------------------

#' Compute the scalar magnitude (scale) adjustment
#'
#' Computes `m = nrow(H) / sum(diag(solve(H, J)))`, the Godambe ratio, where
#' `H` is the negative Hessian (bread) and `J` the centered score outer
#' product (meat) at a single point. If `model` and `data` are supplied a
#' fresh optimizer point is obtained first (with `magnitude_adj = 1.0` and the
#' requested `use_priors`); otherwise the supplied `fit` is used as the point.
#'
#' This is a **single scalar**: it rescales the whole pseudo-likelihood by one
#' global factor (a coarse correction of the overall posterior scale). It is a
#' crude proxy for the full robust covariance and does **not** in general
#' reproduce a direction-dependent multivariate sandwich covariance.
#'
#' @param fit a `CmdStanMCMC` or `CmdStanMLE` object to use as the point when
#'   `model`/`data` are not both supplied.
#' @param use_priors `0L` (MLE) or `1L` (MAP/penalized) -- only used when a
#'   fresh optimizer point is requested via `model`/`data`.
#' @param model an optional `CmdStanModel`; if supplied together with `data`,
#'   a fresh optimizer point is computed.
#' @param data optional data list for a fresh optimizer point.
#' @param init optional starting values (function or list) for the optimizer.
#' @param jacobian logical passed to the Hessian. Defaults to `NULL`, in which
#'   case it is inferred: the optimizer's own Jacobian setting for an MLE
#'   point (matching the reference), `TRUE` otherwise.
#' @param ... further arguments passed to `model$optimize()`.
#' @references
#' Godambe, P. C. (1985). Optimal Robust Inference. Biometrika, 72(3),
#'   611-620.
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
#' # truelik is a TRUE likelihood, so the Godambe ratio is ~ 1:
#' pl_magnitude_adj(fit)
#' @export
pl_magnitude_adj <- function(fit = NULL, use_priors = 0L,
                             model = NULL, data = NULL, init = NULL,
                             jacobian = NULL, ...) {
  if (!is.null(model) && !is.null(data)) {
    opt <- model$optimize(
      data = utils::modifyList(data, list(magnitude_adj = 1.0,
                                          use_priors = use_priors)),
      init = init,
      jacobian = (use_priors == 1L),
      refresh = 0,
      show_messages = FALSE,
      ...
    )
    point <- opt
  } else {
    if (is.null(fit)) {
      stop("provide either `fit`, or both `model` and `data`")
    }
    point <- fit
  }

  if (is.null(jacobian)) {
    # Match the reference get_mag_adj / get_h_unc: for an optimizer point use
    # the Jacobian setting the optimizer itself was run with; otherwise TRUE.
    jacobian <- if (inherits(point, "CmdStanMLE")) {
      isTRUE(point$metadata()$jacobian)
    } else {
      TRUE
    }
  }

  jh <- get_jh(point, j = TRUE, h = TRUE, hinv = FALSE, jacobian = jacobian)
  j <- jh$j
  h <- jh$h
  m <- nrow(h) / sum(diag(solve(h, j)))
  m
}

#' Scalar magnitude (scale) adjustment
#'
#' Thin exported wrapper around [pl_magnitude_adj()] returning the scalar
#' `m` used to rescale a pseudo-likelihood.
#'
#' @inheritParams pl_magnitude_adj
#' @examplesIf identical(Sys.getenv("PP_EXAMPLES_STAN"), "1") && requireNamespace("cmdstanr", quietly = TRUE)
#' # Thin wrapper around [pl_magnitude_adj()]; shown here taking a fresh
#' # optimizer point via model/data rather than a fit:
#' library(cmdstanr)
#' model <- cmdstan_model(
#'   system.file("test-models", "truelik.stan", package = "pseudopost"),
#'   force_recompile = TRUE, quiet = TRUE)
#' set.seed(42)
#' dat <- list(N = 100L, y = rnorm(100, 1, 1),
#'             magnitude_adj = 1.0, use_priors = 1L)
#' magnitude_adjust(model = model, data = dat, init = list(list(theta = 0.0)))
#' @export
magnitude_adjust <- function(fit = NULL, use_priors = 0L,
                             model = NULL, data = NULL, init = NULL,
                             jacobian = NULL, ...) {
  pl_magnitude_adj(fit = fit, use_priors = use_priors, model = model,
                   data = data, init = init, jacobian = jacobian, ...)
}

#' Refit a model with the magnitude adjustment applied
#'
#' Computes the scalar magnitude adjustment, then refits the model with
#' `magnitude_adj = m` and `use_priors = 1L`, returning a no-adjustment
#' `pseudo_post` summary of the refit with the scalar attached as
#' `$magnitude`.
#'
#' @param fit the original `CmdStanMCMC` fit (used for default MCMC settings).
#' @param model a `CmdStanModel` to refit.
#' @param data the data list for the model.
#' @param params character vector of (possibly transformed) parameter names to
#'   summarize; defaults to all base and transformed parameters.
#' @param mcmc_args optional list of arguments for `model$sample()`; defaults
#'   to the chain/iteration/seed settings of `fit`.
#' @param use_priors `0L` (MLE) or `1L` (MAP) for computing the magnitude.
#' @param init optional starting values for the optimizer and the refit.
#' @param ... further arguments passed to `model$sample()`.
#' @references
#' Godambe, P. C. (1985). Optimal Robust Inference. Biometrika, 72(3),
#'   611-620.
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
#' ma <- refit_magnitude_adjust(
#'   fit, model = model, data = dat, init = list(list(theta = 0.0)),
#'   mcmc_args = list(chains = 2L, iter_warmup = 200L, iter_sampling = 200L,
#'                    seed = 42L, refresh = 0, show_messages = FALSE))
#' ma$magnitude
#' as.data.frame(ma)
#' @export
refit_magnitude_adjust <- function(fit, model, data, params = NULL,
                                    mcmc_args = NULL, use_priors = 0L,
                                    init = NULL, ...) {
  # Validate the fit (cmdstanr type + per-observation log-likelihood). The
  # magnitude rescale itself is applied through `model` (the optimizer point
  # and the refit both call `model`), so the MODEL must declare `magnitude_adj`
  # -- not merely the original `fit`. Checking the model catches a fit/model
  # mismatch that a fit-only check would miss.
  check_pl_requirements(fit)
  check_model_data_arg(model, "magnitude_adj")

  m <- pl_magnitude_adj(model = model, data = data,
                        use_priors = use_priors, init = init)

  if (is.null(mcmc_args)) {
    md <- fit$metadata()
    mcmc_args <- list(
      chains = md$num_chains,
      iter_warmup = md$iter_warmup,
      iter_sampling = md$iter_sampling,
      seed = md$seed
    )
  }

  if (is.null(params)) {
    mv <- fit$runset$args$model_variables
    params <- c(names(mv$parameters), names(mv$transformed_parameters))
  }

  f2 <- do.call(
    model$sample,
    c(
      list(data = utils::modifyList(data, list(magnitude_adj = m,
                                               use_priors = 1L)),
           init = init),
      mcmc_args,
      list(...)
    )
  )

  res <- adjust_pseudo_posterior(f2, params, methods = "noadj",
                                 verbose = FALSE)
  res$magnitude <- m
  res
}
