# Initialize cmdstanr model methods where available.
ensure_model_methods <- function(fit) {
  if (inherits(fit, c("CmdStanMCMC", "CmdStanMLE"))) {
    tryCatch(fit$init_model_methods(), error = function(e) NULL)
  }
  invisible(fit)
}

# An MCMC posterior mean is not an MLE or MAP.
resolve_curvadj_point <- function(user_point, model, data, init,
                                  use_priors, jacobian, method) {
  if (!is.null(user_point)) {
    return(user_point)
  }
  if (!is.null(model) && !is.null(data)) {
    return(model$optimize(
      data = utils::modifyList(data, list(magnitude_adj = 1.0,
                                          use_priors = use_priors)),
      init = init,
      jacobian = jacobian,
      refresh = 0,
      show_messages = FALSE
    ))
  }
  which <- if (identical(use_priors, 1L)) "MAP" else "MLE"
  stop(sprintf(
    paste0("%s needs a %s point. Supply `model` and `data` so the %s can be ",
           "optimized (it is called with `use_priors = %d`; for a genuine MLE ",
           "the model must disable ALL priors when `use_priors = 0`), or a ",
           "pre-optimized `opt` (MLE) / `opt_map` (MAP). The MCMC posterior ",
           "mean is not substituted, since it is not a %s."),
    method, which, which, use_priors, which), call. = FALSE)
}

pp_choices <- function(x, choices, name) {
  x <- unname(unlist(x))
  bad <- setdiff(x, choices)
  if (length(bad)) {
    stop(sprintf("%s: unknown value(s) '%s'; must be one of '%s'.",
                 name, paste(bad, collapse = "', '"),
                 paste(choices, collapse = "', '")), call. = FALSE)
  }
  if (!length(x)) x <- choices[1]
  x
}

# Map a sandwich bread to (htype, jacobian, evaluation point).
bread_cfg <- function(bread) {
  switch(bread,
    mcmc        = list(htype = "postvar", jacobian = FALSE, point = "fit"),
    hessian_mle = list(htype = "hessian", jacobian = FALSE, point = "mle"),
    hessian_map = list(htype = "hessian", jacobian = TRUE,  point = "map"))
}

# Systematic label for one (estimator, bread, meat, ci) combination.
combo_label <- function(estimator, bread, meat, ci) {
  cimap <- c(quantile = "quant", wald = "wald")
  if (estimator == "ij") return(paste0("ij.", cimap[ci]))
  bmap <- c(mcmc = "mcmc", hessian_mle = "mle", hessian_map = "map")
  mmap <- c(score = "score", score_avg = "scoreavg")
  pref <- if (estimator == "ofs") "ofs" else "sandwich"
  paste0(pref, ".", bmap[bread], ".", mmap[meat], ".", cimap[ci])
}

# Expand the (estimator, bread, meat, ci) selections into a list of specs.
# `bread`/`meat` are placeholders (ignored) when estimator = "ij".
expand_combos <- function(estimator, bread, meat, ci) {
  est <- pp_choices(estimator, c("ij", "sandwich", "ofs"), "estimator")
  brd <- pp_choices(bread, c("mcmc", "hessian_mle", "hessian_map"), "bread")
  mtt <- pp_choices(meat, c("score", "score_avg"), "meat")
  ci_ <- pp_choices(ci, c("quantile", "wald"), "ci")
  specs <- list()
  add <- function(e, b, m, c) {
    specs[[length(specs) + 1L]] <<- list(
      estimator = e, bread = b, meat = m, ci = c,
      label = combo_label(e, b, m, c))
  }
  for (e in est) {
    if (e == "ij") {
      for (c in ci_) add("ij", "mcmc", "score", c)
    } else {
      for (b in brd) for (m in mtt) for (c in ci_) add(e, b, m, c)
    }
  }
  specs
}

#' Adjust a pseudo-posterior's uncertainty
#'
#' Runs one or more post-hoc covariance (sandwich) adjustments on a fit to a
#' composite / pairwise / pseudo-likelihood Stan model and returns corrected
#' standard errors, credible intervals, and adjusted posterior draws.
#'
#' @details
#' The robust covariance `V = H^-1 J H^-1` (Shaby, 2014) is estimated by one of
#' three estimators, selected by `estimator`:
#' \describe{
#'   \item{`"ij"`}{The **infinitesimal jackknife (IJ)** adjustment: a
#'     model-blind estimate of the full `V` straight from the MCMC draws
#'     ([pl_ij()]). `bread`/`meat` are not used.}
#'   \item{`"sandwich"`}{The **closed** Godambe sandwich `H^-1 J H^-1`, with an
#'     explicit precision `bread` (`H`) and score `meat` (`J`).}
#'   \item{`"ofs"`}{The **open-faced** sandwich (Shaby, 2014): computes
#'     `Omega = (Q^-1 P)^(1/2)` from an explicit precision `bread` (`Q`) and
#'     score `meat` (`P`) and applies it directly to the centered draws
#'     without calibrating to the empirical draw covariance. Assumes the
#'     draws' covariance is the bread `Q^-1` (the "one-bread" asymptotic).
#'     Uses `bread`/`meat` like `"sandwich"`.}
#' }
#' Under `estimator = "sandwich"` or `estimator = "ofs"`:
#' \describe{
#'   \item{`bread`}{`"mcmc"` (inverse MCMC draw covariance), `"hessian_mle"`
#'     (negative Hessian of the pseudo-likelihood at the MLE, no Jacobian), or
#'     `"hessian_map"` (negative Hessian of the pseudo-posterior at the MAP,
#'     with Jacobian).}
#'   \item{`meat`}{`"score"` (the moment estimator from the
#'     per-observation score vector) or `"score_avg"` (a draw-averaged version).}
#' }
#' `ci` selects the interval for every estimator: `"quantile"` (percentile of
#' the adjusted draws) or `"wald"` (estimate +/- z * se).
#'
#' The `"hessian_mle"`/`"hessian_map"` breads evaluate the Hessian at a
#' (penalized-)likelihood optimum, so a `model`/`data` pair or an explicit
#' `opt`/`opt_map` is required; there is no MCMC-posterior-mean fallback,
#' since the posterior mean is not an MLE or MAP and would mislabel the point.
#' For `"hessian_mle"` to be a true likelihood-only MLE, the model must
#' implement a mode that disables all priors (e.g. via `use_priors = 0`); the
#' package cannot verify that the Stan program honors that switch.
#'
#' All of `estimator`, `bread`, `meat`, and `ci` may be vectors to run the
#' corresponding sub-grid; each valid combination produces one row per
#' requested parameter.
#'
#' **Legacy interface:** passing `methods` selects the original named presets,
#' which are silent aliases for the axes above: `"noadj"` (baseline),
#' `"ijse"`/`"ijse2"` (IJ, Wald/quantile), `"curvadj"`/`"curvadj2"`/
#' `"curvadj3"` (sandwich with `bread` = `hessian_mle`/`hessian_map`/`mcmc`,
#' `meat` = `score`, quantile). When `methods` is supplied it takes precedence
#' over `estimator`/`bread`/`meat`/`ci`, and `ij_method`/`jacobian` apply.
#'
#' @param fit a `CmdStanMCMC` (or `CmdStanMLE`) object that satisfies the
#'   pseudo-likelihood requirements (see [check_pl_requirements()]).
#' @param params character vector of (possibly transformed) parameter names to
#'   report on.
#' @param methods optional character vector of legacy presets (see Details).
#'   `NULL` (the default) uses the `estimator`/`bread`/`meat`/`ci` interface.
#' @param estimator `"ij"` (infinitesimal jackknife, default), `"sandwich"`
#'   (closed Godambe), or `"ofs"` (open-faced sandwich). A vector runs the
#'   corresponding sub-grid.
#' @param bread the precision `H` for `"sandwich"` and `"ofs"`: `"mcmc"`,
#'   `"hessian_mle"`, or `"hessian_map"`. `"hessian_mle"` evaluates the
#'   Hessian at a likelihood-only MLE (the model must disable ALL priors when
#'   `use_priors = 0`) and requires `model`/`data` (or a pre-optimized
#'   `opt`); `"hessian_map"` evaluates it at the MAP and requires
#'   `model`/`data` (or `opt_map`). Neither falls back to the MCMC posterior
#'   mean. Ignored for `"ij"`.
#' @param meat the `J` for `"sandwich"` and `"ofs"`: `"score"` or
#'   `"score_avg"`. Ignored for `"ij"`.
#' @param ci `"quantile"` (default) or `"wald"`; applies to both estimators.
#' @param name_lli generated-quantity name for per-observation
#'   log-likelihood contributions.
#' @param ij_method legacy: the IJ variant, `"proxy"` (default) or `"score"`;
#'   used only when `methods` selects `ijse`/`ijse2`.
#' @param score_draws number of draws subsampled for score-based meat.
#' @param ci_level credible-interval level (default `0.90`).
#' @param model optional `CmdStanModel`; with `data`, used to (re)optimize the
#'   sandwich evaluation points.
#' @param data optional data list, paired with `model`.
#' @param opt a pre-optimized likelihood-only MLE point (`CmdStanMLE`) used
#'   for `bread = "hessian_mle"`; it must be a point at which all priors are
#'   disabled.
#' @param opt_map a pre-optimized MAP point (`CmdStanMLE`) used for
#'   `bread = "hessian_map"`.
#' @param init optional starting values for any optimizer calls.
#' @param jacobian legacy; reserved for API stability (the bread's Jacobian
#'   setting is now fixed by `bread`).
#' @param verbose if `TRUE`, print a message when a method fails.
#' @return an object of class `pseudo_post`: a list with
#'   \itemize{
#'     \item `table` -- a data frame (long form) with columns
#'       `method, param, est, se, lo, hi`.
#'     \item `draws` -- a named list of adjusted-draws matrices, one per method
#'       that produced draws.
#'     \item `meta` -- a list with `fit`, `methods`, `ci_level`, `ij_method`,
#'       `warnings` and `notes`.
#'   }
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
#' pp <- adjust_pseudo_posterior(
#'   fit, "theta",
#'   methods = c("noadj", "ijse", "ijse2", "curvadj", "curvadj2", "curvadj3"),
#'   model = model, data = dat, init = list(list(theta = 0.0)),
#'   verbose = FALSE)
#' print(pp)
#' @export
adjust_pseudo_posterior <- function(fit, params,
                                    methods = NULL,
                                    estimator = "ij",
                                    bread = "mcmc",
                                    meat = "score",
                                    ci = "quantile",
                                    name_lli = "log_lik",
                                    ij_method = "proxy",
                                    score_draws = 200L,
                                    ci_level = 0.90,
                                    model = NULL, data = NULL,
                                    opt = NULL, opt_map = NULL,
                                    init = NULL, jacobian = TRUE,
                                    verbose = TRUE) {
  check_pl_requirements(fit, name_lli = name_lli)

  alpha <- (1 - ci_level) / 2
  lo_q <- alpha
  hi_q <- 1 - alpha
  z <- qnorm(1 - alpha)

  if (is.null(methods) || any(methods != "noadj")) {
    ensure_model_methods(fit)
  }

  rows <- list()
  draws_out <- list()
  warnings_out <- character(0)
  notes <- character(0)

  quantile_ci <- function(d) {
    list(est = colMeans(d),
         se = apply(d, 2, sd),
         lo = apply(d, 2, quantile, probs = lo_q),
         hi = apply(d, 2, quantile, probs = hi_q))
  }

  result_from_draws <- function(d) {
    ci <- quantile_ci(d)
    list(draws = d, est = ci$est, se = ci$se, lo = ci$lo, hi = ci$hi)
  }

  ci_res <- function(d, ci) {
    if (ci == "wald") {
      est <- colMeans(d)
      se <- apply(d, 2, sd)
      list(draws = d, est = est, se = se,
           lo = est - z * se, hi = est + z * se)
    } else {
      result_from_draws(d)
    }
  }

  if (!is.null(methods)) {
    methods <- match.arg(
      methods, c("noadj", "ijse", "ijse2", "curvadj", "curvadj2", "curvadj3"),
      several.ok = TRUE)
    ij_method <- match.arg(ij_method, c("proxy", "score"))
    requested <- methods

    for (m in methods) {
      res <- tryCatch(
        switch(m,
          noadj = {
            d <- fit$draws(params, format = "draws_matrix")
            result_from_draws(d)
          },
          ijse = {
            B <- pl_ij(fit, name_lli = name_lli, method = ij_method,
                              score_draws = score_draws)
            d <- transform_to_target(fit, params, T = B, name_lli = name_lli)
            est <- colMeans(d)
            se <- apply(d, 2, sd)
            list(draws = d, est = est, se = se,
                 lo = est - z * se, hi = est + z * se)
          },
          ijse2 = {
            B <- pl_ij(fit, name_lli = name_lli, method = ij_method,
                              score_draws = score_draws)
            d <- transform_to_target(fit, params, T = B, name_lli = name_lli)
            result_from_draws(d)
          },
          curvadj = {
            point <- resolve_curvadj_point(opt, model, data, init,
                                           use_priors = 0L, jacobian = FALSE,
                                           method = "curvadj")
            g <- pl_godambe(point, htype = "hessian", name_lli = name_lli,
                            jacobian = FALSE)
            d <- transform_to_target(fit, params, T = g$V, name_lli = name_lli)
            result_from_draws(d)
          },
          curvadj2 = {
            point <- resolve_curvadj_point(opt_map, model, data, init,
                                           use_priors = 1L, jacobian = TRUE,
                                           method = "curvadj2")
            g <- pl_godambe(point, htype = "hessian", name_lli = name_lli,
                            jacobian = TRUE)
            d <- transform_to_target(fit, params, T = g$V, name_lli = name_lli)
            result_from_draws(d)
          },
          curvadj3 = {
            g <- pl_godambe(fit, htype = "postvar", name_lli = name_lli)
            d <- transform_to_target(fit, params, T = g$V, name_lli = name_lli)
            result_from_draws(d)
          }
        ),
        error = function(e) {
          warnings_out <<- c(warnings_out,
                             paste0(m, ": ", conditionMessage(e)))
          if (verbose) {
            message("method '", m, "' failed: ", conditionMessage(e))
          }
          NULL
        }
      )

      if (!is.null(res)) {
        dmat <- res$draws
        rows[[m]] <- data.frame(
          method = m,
          param = colnames(dmat),
          est = unname(res$est),
          se = unname(res$se),
          lo = unname(res$lo),
          hi = unname(res$hi),
          stringsAsFactors = FALSE
        )
        draws_out[[m]] <- dmat
      }
    }
  } else {
    specs <- expand_combos(estimator, bread, meat, ci)
    requested <- vapply(specs, function(s) s$label, character(1))
    draws_cache <- new.env(parent = emptyenv())
    point_cache <- new.env(parent = emptyenv())
    avg_meat <- NULL

    for (spec in specs) {
      lab <- spec$label
      res <- tryCatch({
        key <- if (spec$estimator == "ij") "ij" else
          paste0(spec$estimator, ".", spec$bread, ".", spec$meat)
        if (exists(key, envir = draws_cache, inherits = FALSE)) {
          d <- get(key, envir = draws_cache, inherits = FALSE)
        } else {
          if (spec$estimator == "ij") {
            T <- pl_ij(fit, name_lli = name_lli, method = "proxy",
                               score_draws = score_draws)
            d <- transform_to_target(fit, params, T = T, name_lli = name_lli)
          } else if (spec$estimator == "ofs") {
            cfg <- bread_cfg(spec$bread)
            jmeat <- NULL
            if (spec$meat == "score_avg") {
              if (is.null(avg_meat)) {
                avg_meat <- get_j_avg(fit, score_draws = score_draws,
                                      name_lli = name_lli)
              }
              jmeat <- avg_meat
            }
            if (cfg$point == "fit") {
              point <- fit
            } else {
              if (!exists(cfg$point, envir = point_cache, inherits = FALSE)) {
                user_point <- if (cfg$point == "map") opt_map else opt
                point <- resolve_curvadj_point(
                  user_point, model, data, init,
                  use_priors = if (cfg$point == "map") 1L else 0L,
                  jacobian = cfg$jacobian, method = lab)
                assign(cfg$point, point, envir = point_cache)
              }
              point <- get(cfg$point, envir = point_cache, inherits = FALSE)
            }
            g <- pl_ofs(point, htype = cfg$htype, name_lli = name_lli,
                        jacobian = cfg$jacobian, j = jmeat)
            d <- transform_ofs(fit, params, Omega = g$Omega,
                               name_lli = name_lli)
          } else {
            cfg <- bread_cfg(spec$bread)
            jmeat <- NULL
            if (spec$meat == "score_avg") {
              if (is.null(avg_meat)) {
                avg_meat <- get_j_avg(fit, score_draws = score_draws,
                                      name_lli = name_lli)
              }
              jmeat <- avg_meat
            }
            if (cfg$point == "fit") {
              point <- fit
            } else {
              if (!exists(cfg$point, envir = point_cache, inherits = FALSE)) {
                user_point <- if (cfg$point == "map") opt_map else opt
                point <- resolve_curvadj_point(
                  user_point, model, data, init,
                  use_priors = if (cfg$point == "map") 1L else 0L,
                  jacobian = cfg$jacobian, method = lab)
                assign(cfg$point, point, envir = point_cache)
              }
              point <- get(cfg$point, envir = point_cache, inherits = FALSE)
            }
            g <- pl_godambe(point, htype = cfg$htype, name_lli = name_lli,
                            jacobian = cfg$jacobian, j = jmeat)
            T <- g$V
            d <- transform_to_target(fit, params, T = T, name_lli = name_lli)
          }
          assign(key, d, envir = draws_cache)
        }
        ci_res(d, spec$ci)
      },
      error = function(e) {
        warnings_out <<- c(warnings_out,
                           paste0(lab, ": ", conditionMessage(e)))
        if (verbose) {
          message("method '", lab, "' failed: ", conditionMessage(e))
        }
        NULL
      }
      )

      if (!is.null(res)) {
        dmat <- res$draws
        rows[[lab]] <- data.frame(
          method = lab,
          param = colnames(dmat),
          est = unname(res$est),
          se = unname(res$se),
          lo = unname(res$lo),
          hi = unname(res$hi),
          stringsAsFactors = FALSE
        )
        draws_out[[lab]] <- dmat
      }
    }
  }

  table <- do.call(rbind, rows)
  if (is.null(table)) {
    table <- data.frame(method = character(0), param = character(0),
                        est = numeric(0), se = numeric(0),
                        lo = numeric(0), hi = numeric(0),
                        stringsAsFactors = FALSE)
  }

  structure(
    list(table = table, draws = draws_out,
         meta = list(fit = fit, methods = requested, ci_level = ci_level,
                     ij_method = ij_method, warnings = warnings_out,
                     notes = notes)),
    class = "pseudo_post"
  )
}
