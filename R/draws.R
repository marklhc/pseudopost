#' Affinely adjust posterior draws to a target covariance
#'
#' Centers the draws at their column means, post-multiplies by `t(M)`, and
#' re-centers. If the draws have covariance `C`, the adjusted draws have
#' covariance `M C M'`. Setting `M = mat_sqrt(T) %*% mat_inv_sqrt(C)` maps
#' `C` to `T` (exact when `C` is the empirical draw covariance).
#'
#' @param draws_mat a matrix of draws (rows = iterations, columns = parameters
#'   in the unconstrained space).
#' @param M the transformation matrix.
#' @return a matrix of adjusted draws with the same dimensions and column
#'   names as `draws_mat`.
#' @examples
#' set.seed(1)
#' d <- matrix(rnorm(2000 * 2), 2000, 2)
#' colnames(d) <- c("a", "b")
#' M <- diag(c(2, 1))             # double the spread of "a"
#' d_adj <- adjust_posterior_draws(d, M)
#' round(cov(d), 2)               # ~ identity
#' round(cov(d_adj), 2)           # ~ diag(4, 1)
#' @export
adjust_posterior_draws <- function(draws_mat, M) {
  theta_hat <- colMeans(draws_mat)
  draws_centered <- sweep(draws_mat, 2, theta_hat, "-")
  draws_adj_centered <- draws_centered %*% t(M)
  draws_adj <- sweep(draws_adj_centered, 2, theta_hat, "+")
  colnames(draws_adj) <- colnames(draws_mat)
  draws_adj
}

#' Extract a (possibly bracketed) scalar from a constrained vector
#'
#' `fit$constrain_variables()` returns a list keyed by TOP-LEVEL names only
#' (e.g. `"lambda_delta"` -> a vector, `"thres_delta"` -> an array), whereas
#' `fit$draws(params)` reports each scalar column with a bracketed name
#' (e.g. `"thres_delta[2,1]"`). This helper resolves a requested scalar column
#' by parsing its base name and integer subscripts, rather than by direct list
#' lookup.
#'
#' @param cc a list of top-level constrained values as returned by
#'   `fit$constrain_variables()`.
#' @param name a (possibly bracketed) parameter name, e.g. `"thres_delta[2,1]"`.
#' @return a numeric scalar.
#' @keywords internal
extract_constrained_scalar <- function(cc, name) {
  if (grepl("[", name, fixed = TRUE)) {
    base  <- sub("\\[.*$", "", name)
    idx   <- sub("^.*\\[", "", name)
    idx   <- sub("\\]$", "", idx)
    parts <- as.integer(strsplit(idx, ",")[[1]])
    # Explicit per-dimension indexing: m[1, 1] (NOT m[c(1, 1)], which R reads
    # as row indices in column 1 and returns a length-2 vector).
    do.call("[", c(list(cc[[base]]), as.list(parts)))
  } else {
    cc[[name]]
  }
}

#' Constrain a matrix of unconstrained draws row-by-row
#'
#' Each row of `unconstrained_draws` is passed through
#' `fit$constrain_variables()` to collect the requested parameters.
#'
#' @param fit a `CmdStanMCMC` or `CmdStanMLE` object.
#' @param unconstrained_draws a matrix of unconstrained parameter draws
#'   (rows = iterations).
#' @param params character vector of (possibly transformed) parameter names to
#'   extract from the constrained vector.
#' @return a numeric matrix of constrained draws (rows = iterations).
#' @keywords internal
constrain_draws <- function(fit, unconstrained_draws, params) {
  out <- fit$draws(params, format = "draws_matrix")
  out[] <- NA_real_
  cols <- colnames(out)
  for (i in seq_len(nrow(unconstrained_draws))) {
    cc <- fit$constrain_variables(unconstrained_draws[i, ])
    for (j in seq_along(cols)) {
      out[i, j] <- extract_constrained_scalar(cc, cols[j])
    }
  }
  out
}

#' Transform posterior draws so their covariance matches a target
#'
#' The unifying transform behind every adjustment family. The unconstrained
#' draws have empirical covariance `C = cov(fit$unconstrain_draws())`; the map
#' `M = mat_sqrt(T) %*% mat_inv_sqrt(C)` is applied to the centered draws
#' (`adjust_posterior_draws`) so their covariance in the unconstrained space
#' becomes exactly `T`, then the draws are mapped back to the requested
#' constrained parameters (`constrain_draws`).
#'
#' The map is built from the ACTUAL draw covariance `C`, so the adjusted draws
#' have covariance `T` in the unconstrained space for ANY target `T` (this is
#' what makes the reported SEs/quantiles consistent with the target, rather
#' than with a point's inverse-Hessian, which may differ from the draws'
#' empirical covariance). The exactness holds in the unconstrained space;
#' mapping back to the constrained `params` can change the covariance in the
#' constrained space (the affine map is not preserved under non-linear
#' transforms).
#'
#' @param fit a `CmdStanMCMC` (or `CmdStanMLE`) object (source of the draws).
#' @param params character vector of (possibly transformed) parameter names to
#'   return.
#' @param T the target covariance matrix in the unconstrained space.
#' @param name_lli generated-quantity name for per-observation
#'   log-likelihood contributions (unused by the transform itself; accepted
#'   for signature consistency with the adjustment entry point).
#' @return a numeric matrix of adjusted, constrained draws
#'   (rows = iterations, columns = `params`).
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
#' u <- fit$unconstrain_draws(format = "draws_matrix")
#' C_emp <- cov(u)
#' # T = C_emp makes the transform the identity (up to round-trip noise):
#' out <- transform_to_target(fit, "theta", T = C_emp)
#' head(out)
#' @export
transform_to_target <- function(fit, params, T, name_lli = "log_lik") {
  u <- fit$unconstrain_draws(format = "draws_matrix")
  C_emp <- cov(u)
  M <- mat_sqrt(T) %*% mat_inv_sqrt(C_emp)
  u_adj <- adjust_posterior_draws(u, M)
  constrain_draws(fit, u_adj, params)
}

#' Apply the OFS draw transform (no C_emp calibration)
#'
#' The open-faced sandwich variant of [transform_to_target()]: applies the
#' square root `Omega` from [pl_ofs()] directly to the centered unconstrained
#' draws (`adjust_posterior_draws` with `M = Omega`) and maps back to the
#' requested constrained parameters, without calibrating to the empirical draw
#' covariance `C_emp`.
#'
#' @param fit a `CmdStanMCMC` (or `CmdStanMLE`) object (source of the draws).
#' @param params character vector of (possibly transformed) parameter names to
#'   return.
#' @param Omega the open-faced sandwich square root from [pl_ofs()] (the
#'   transformation matrix applied to the centered draws).
#' @param name_lli generated-quantity name for per-observation
#'   log-likelihood contributions (unused by the transform itself; accepted
#'   for signature consistency with the adjustment entry point).
#' @return a numeric matrix of adjusted, constrained draws
#'   (rows = iterations, columns = `params`).
#' @keywords internal
transform_ofs <- function(fit, params, Omega, name_lli = "log_lik") {
  u <- fit$unconstrain_draws(format = "draws_matrix")
  u_adj <- adjust_posterior_draws(u, Omega)
  constrain_draws(fit, u_adj, params)
}
