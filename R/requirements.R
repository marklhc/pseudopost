# ---------------------------------------------------------------------------
# requirements.R
#
# Pseudo-likelihood requirements validation. A model is eligible for post-hoc
# adjustment if (a) it is a cmdstanr fit, and (b) it exposes a
# per-observation log-likelihood generated quantity (a vector, one entry per
# pseudo-observation). When magnitude adjustment is requested the model must
# additionally take a `magnitude_adj` data argument.
# ---------------------------------------------------------------------------

#' Check the pseudo-likelihood requirements of a fit
#'
#' Verifies that a cmdstanr fit exposes a per-observation log-likelihood
#' generated quantity (required for any adjustment) and, optionally, a
#' `magnitude_adj` data argument (required for magnitude adjustment).
#'
#' @details
#' This check validates only that `log_lik` exists and is vector-shaped (one
#' entry per pseudo-observation); it does NOT verify that its entries are
#' contributions to independent observational units or that they sum to the
#' fitted pseudo-likelihood -- those remain the caller's (model author's)
#' responsibility. A scalar `log_lik` (a single generated quantity, i.e. one
#' pseudo-observation) is accepted with a warning that the IJ/sandwich
#' correction is then degenerate, rather than being rejected as
#' non-vector-shaped.
#'
#' @param fit a `CmdStanMCMC` or `CmdStanMLE` object.
#' @param name_lli name of the generated quantity holding per-observation
#'   log-likelihood contributions (default `"log_lik"`).
#' @param need_magnitude if `TRUE`, also require a `magnitude_adj` data
#'   argument in the model.
#' @return invisibly, a list `list(m, ok = TRUE)` where `m` is the length of
#'   the `log_lik` vector (the number of pseudo-observations). Signals an
#'   error if the requirements are not met.
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
#' check_pl_requirements(fit)                 # m = 100
#' check_pl_requirements(fit, need_magnitude = TRUE)
#' @export
check_pl_requirements <- function(fit, name_lli = "log_lik", need_magnitude = FALSE) {
  # The adjustment machinery drives the fit through cmdstanr's model-method
  # interface (hessian / constrain_variables / ...); make the runtime
  # dependency explicit here at the validation boundary.
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("pseudopost requires the 'cmdstanr' package to inspect Stan fits; ",
         "install it with install.packages('cmdstanr').", call. = FALSE)
  }
  if (!inherits(fit, c("CmdStanMCMC", "CmdStanMLE"))) {
    stop("`fit` must be a CmdStanMCMC or CmdStanMLE object, not '",
         paste(class(fit), collapse = "/"), "'.")
  }

  lli <- tryCatch(
    fit$draws(name_lli, format = "draws_matrix"),
    error = function(e) {
      stop("generated quantity '", name_lli, "' not found in fit; the model ",
           "must expose a per-observation log-likelihood vector in ",
            "`generated quantities` (the pseudo-likelihood requirements).",
           call. = FALSE)
    }
  )
  # The requirements call for a VECTOR (one entry per pseudo-observation). A
  # matrix-valued generated quantity is stored as scalar names with 2-D
  # brackets (e.g. `log_lik[i, j]`); its flattened length would be mistaken
  # for the number of pseudo-observations, so reject it explicitly. (Draw
  # matrices with no column names -- e.g. from minimal mock fits -- are left
  # to the length check.)
  cn <- colnames(lli)
  if (!is.null(cn) && any(grepl("\\[[^\\]]*,", cn))) {
    stop("generated quantity '", name_lli, "' appears to be matrix-valued ",
          "(e.g. '", name_lli, "[i, j]'). The pseudo-likelihood requirements",
          " call for a VECTOR with one entry per pseudo-observation; flatten ",
         "the per-observation contributions into a vector.", call. = FALSE)
  }
  m <- ncol(lli)

  if (m == 1) {
    warning("only 1 pseudo-observation; IJ/sandwich correction is ",
            "degenerate", call. = FALSE)
  }

  if (need_magnitude) {
    # Where the declared data arguments live depends on the cmdstanr version:
    # newer versions may expose them via metadata(); in 0.9.0 they are in the
    # parsed model variables stored on the runset.
    md <- fit$metadata()
    has_mag <- FALSE
    if (!is.null(md$data)) {
      has_mag <- "magnitude_adj" %in% names(md$data)
    }
    if (!has_mag && !is.null(fit$runset)) {
      mv <- fit$runset$args$model_variables
      if (!is.null(mv) && !is.null(mv$data)) {
        has_mag <- "magnitude_adj" %in% names(mv$data)
      }
    }
    if (!has_mag) {
      stop("MagAdj requires a 'magnitude_adj' data argument in the model",
           call. = FALSE)
    }
  }

  invisible(list(m = m, ok = TRUE))
}

#' Check that a Stan model declares a named data argument
#'
#' Validates that a `CmdStanModel` declares the given data argument by
#' inspecting the model's declared data variables (`model$variables()$data`).
#' This is the model-side counterpart of the `need_magnitude` check in
#' [check_pl_requirements()]: it operates on the `CmdStanModel` that a quantity
#' will actually be supplied to (e.g. the refit model), rather than on an
#' existing fit, so a mismatch between a fit and its refit model is caught.
#'
#' @param model a `CmdStanModel` object.
#' @param name the data-argument name to require (default `"magnitude_adj"`).
#' @return invisibly `TRUE`. Signals an error if `model` does not declare
#'   `name`.
#' @keywords internal
check_model_data_arg <- function(model, name = "magnitude_adj") {
  declared <- tryCatch(names(model$variables()$data), error = function(e) NULL)
  if (is.null(declared)) {
    stop("could not read the model's declared data arguments; `model` must ",
         "be a compiled cmdstanr CmdStanModel.", call. = FALSE)
  }
  if (!name %in% declared) {
    stop("the model does not declare a '", name, "' data argument, which is ",
         "required to apply the magnitude adjustment.", call. = FALSE)
  }
  invisible(TRUE)
}
