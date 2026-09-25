# ---------------------------------------------------------------------------
# Requirements tests for check_pl_requirements().
# The pure-R parts (mock fits) always run; the real-fit parts are gated by
# skip_stan() (see helper-fit.R).
# ---------------------------------------------------------------------------

# A minimal CmdStanMCMC stand-in: an environment with the two accessors
# check_pl_requirements uses (draws, metadata; plus an optional runset for the
# cmdstanr 0.9.0 data-argument lookup path).
make_mock_fit <- function(gq = NULL, metadata = NULL,
                          has_runset = FALSE, runset_data = NULL) {
  e <- new.env(parent = emptyenv())
  e$draws <- function(name, format = "draws_matrix", ...) {
    if (is.null(gq) || is.null(gq[[name]])) {
      stop("no such variable: ", name)
    }
    gq[[name]]
  }
  e$metadata <- function() metadata
  if (has_runset) {
    e$runset <- list(args = list(model_variables = list(data = runset_data)))
  }
  class(e) <- "CmdStanMCMC"
  e
}

test_that("check_pl_requirements rejects non-cmdstanr objects", {
  expect_error(check_pl_requirements(list()), "CmdStanMCMC or CmdStanMLE")
  expect_error(check_pl_requirements("not a fit"), "CmdStanMCMC or CmdStanMLE")
})

test_that("check_pl_requirements errors when the generated quantity is missing", {
  # Pure-R: a mock whose draws() cannot find the requested quantity.
  mock <- make_mock_fit(gq = list(other_gq = matrix(0, 100, 40)))
  expect_error(check_pl_requirements(mock, name_lli = "does_not_exist"),
               "generated quantity 'does_not_exist' not found")
})

test_that("check_pl_requirements passes and reports m for a valid mock fit", {
  # 100 draws x 40 pseudo-observations.
  mock <- make_mock_fit(gq = list(log_lik = matrix(0, 100, 40)))
  res <- check_pl_requirements(mock)
  expect_true(res$ok)
  expect_equal(res$m, 40)
})

test_that("check_pl_requirements warns when there is a single pseudo-observation", {
  mock <- make_mock_fit(gq = list(log_lik = matrix(0, 100, 1)))
  expect_warning(res <- check_pl_requirements(mock), "degenerate")
  expect_equal(res$m, 1)
})

test_that("check_pl_requirements(need_magnitude = TRUE) requires magnitude_adj", {
  lli <- matrix(0, 100, 40)
  # Declared via metadata$data:
  ok_meta <- make_mock_fit(gq = list(log_lik = lli),
                           metadata = list(data = list(magnitude_adj = 1.0)))
  expect_no_error(check_pl_requirements(ok_meta, need_magnitude = TRUE))
  # Declared via runset$args$model_variables$data (cmdstanr 0.9.0 path):
  ok_runset <- make_mock_fit(gq = list(log_lik = lli),
                             metadata = list(data = NULL),
                             has_runset = TRUE,
                             runset_data = list(magnitude_adj = 1.0, N = 40L))
  expect_no_error(check_pl_requirements(ok_runset, need_magnitude = TRUE))
  # Declared in neither:
  bad <- make_mock_fit(gq = list(log_lik = lli), metadata = list(data = NULL))
  expect_error(check_pl_requirements(bad, need_magnitude = TRUE), "magnitude_adj")
})

# --- Stan-backed requirements checks (gated) --------------------------------

test_that("check_pl_requirements works on a real truelik fit", {
  skip_stan()
  fit <- get_truelik_fit()
  res <- check_pl_requirements(fit)
  expect_true(res$ok)
  expect_equal(res$m, 100) # vector[N] log_lik, N = 100
  expect_error(check_pl_requirements(fit, name_lli = "does_not_exist"),
               "generated quantity 'does_not_exist' not found")
})

test_that("check_pl_requirements(need_magnitude = TRUE) passes for tinypair", {
  skip_stan()
  fit <- get_tinypair_fit()
  res <- check_pl_requirements(fit, need_magnitude = TRUE)
  expect_true(res$ok)
  expect_equal(res$m, 40) # vector[N] log_lik, N = 40
})

# --- Vector shape checks with named draws columns (always run) --------------
# The existing mocks above return raw matrices with NULL colnames, so the
# bracketed-name shape check is skipped for them. These mocks add column
# names to exercise that check directly.

test_that("check_pl_requirements accepts a 1-D bracketed vector log_lik[i] and reports m", {
  lli <- matrix(0, 100, 40)
  colnames(lli) <- paste0("log_lik[", 1:40, "]")
  mock <- make_mock_fit(gq = list(log_lik = lli))
  res <- check_pl_requirements(mock)
  expect_true(res$ok)
  expect_equal(res$m, 40)
})

test_that("check_pl_requirements rejects a matrix-valued log_lik[i, j]", {
  lli <- matrix(0, 100, 4)
  colnames(lli) <- c("log_lik[1,1]", "log_lik[1,2]",
                     "log_lik[2,1]", "log_lik[2,2]")
  mock <- make_mock_fit(gq = list(log_lik = lli))
  expect_error(check_pl_requirements(mock), "matrix-valued")
})

