# ---------------------------------------------------------------------------
# s3.R
#
# S3 methods for the `pseudo_post` class.
# ---------------------------------------------------------------------------

# round() only the numeric columns of a data frame (base round.data.frame
# chokes on the character method/param columns).
round_df <- function(df, digits = 4) {
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], round, digits = digits)
  df
}

#' Print a `pseudo_post` object
#'
#' Prints a tidy summary (rounded) of the adjustment table plus the methods,
#' parameters, CI level, and any warnings/notes.
#'
#' @param x a `pseudo_post` object.
#' @param ... ignored.
#' @examples
#' pp <- structure(
#'   list(table = data.frame(method = "noadj", param = "theta",
#'                           est = 1.03, se = 0.09, lo = 0.86, hi = 1.20,
#'                           stringsAsFactors = FALSE),
#'        draws = list(noadj = matrix(rnorm(20), 10, 2)),
#'        meta = list(fit = NULL, methods = "noadj", ci_level = 0.90,
#'                    ij_method = "proxy", warnings = character(0),
#'                    notes = character(0))),
#'   class = "pseudo_post")
#' print(pp)
#' @export
print.pseudo_post <- function(x, ...) {
  cat("Pseudo-posterior uncertainty adjustment\n")
  cat("Methods   :", paste(x$meta$methods, collapse = ", "), "\n")
  cat("CI level  :", x$meta$ci_level, "\n")
  cat("IJ method:", x$meta$ij_method, "\n")
  cat("Params    :", paste(unique(x$table$param), collapse = ", "), "\n\n")

  tab <- x$table
  if (nrow(tab) > 0) {
    print(round_df(tab, digits = 4), row.names = FALSE)
  } else {
    cat("(no methods succeeded)\n")
  }

  msgs <- c(x$meta$warnings, x$meta$notes)
  if (length(msgs) > 0) {
    cat("\nWarnings / notes:\n")
    for (msg in msgs) cat(" - ", msg, "\n", sep = "")
  }
  invisible(x)
}

#' Summarize a `pseudo_post` object
#'
#' Returns (and prints) the full summary table with a header and any
#' per-method warnings/notes.
#'
#' @param object a `pseudo_post` object.
#' @param ... ignored.
#' @examples
#' pp <- structure(
#'   list(table = data.frame(method = "noadj", param = "theta",
#'                           est = 1.03, se = 0.09, lo = 0.86, hi = 1.20,
#'                           stringsAsFactors = FALSE),
#'        draws = list(noadj = matrix(rnorm(20), 10, 2)),
#'        meta = list(fit = NULL, methods = "noadj", ci_level = 0.90,
#'                    ij_method = "proxy", warnings = character(0),
#'                    notes = character(0))),
#'   class = "pseudo_post")
#' summary(pp)
#' @export
summary.pseudo_post <- function(object, ...) {
  cat("== pseudo_post summary ==\n")
  cat("Methods  :", paste(object$meta$methods, collapse = ", "), "\n")
  cat("CI level :", object$meta$ci_level,
      " (alpha =", round((1 - object$meta$ci_level) / 2, 4), ")\n")
  msgs <- c(object$meta$warnings, object$meta$notes)
  if (length(msgs) > 0) {
    cat("Notes:\n")
    for (msg in msgs) cat("  - ", msg, "\n", sep = "")
  }
  cat("\n")
  print(round_df(object$table, digits = 4), row.names = FALSE)
  invisible(object$table)
}

#' Coerce a `pseudo_post` object to a data frame
#'
#' @param x a `pseudo_post` object.
#' @param ... ignored.
#' @examples
#' pp <- structure(
#'   list(table = data.frame(method = "noadj", param = "theta",
#'                           est = 1.03, se = 0.09, lo = 0.86, hi = 1.20,
#'                           stringsAsFactors = FALSE),
#'        draws = list(noadj = matrix(rnorm(20), 10, 2)),
#'        meta = list(fit = NULL, methods = "noadj", ci_level = 0.90,
#'                    ij_method = "proxy", warnings = character(0),
#'                    notes = character(0))),
#'   class = "pseudo_post")
#' as.data.frame(pp)
#' @export
as.data.frame.pseudo_post <- function(x, ...) {
  x$table
}
