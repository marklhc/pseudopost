# ---------------------------------------------------------------------------
# Tests for the scalar magnitude (scale) adjustment entry points.
#
# The magnitude rescale is applied through the REFIT model, so that model (not
# merely the original fit) must declare a `magnitude_adj` data argument. These
# tests pin down that a model WITHOUT that argument is rejected up front
# (before the optimize/refit path), including the fit/model-mismatch case where
# the fit has the argument but the refit model does not. Stan-backed (gated by
# skip_stan(); see helper-fit.R).
# ---------------------------------------------------------------------------

test_that("refit_magnitude_adjust rejects a model without a magnitude_adj argument", {
  skip_stan()
  fit <- get_nomag_fit()
  # The requirement check fires before any optimization/refit, so this errors
  # quickly and with the specific magnitude_adj message (not a Stan error).
  expect_error(
    refit_magnitude_adjust(fit, model = get_nomag_model(),
                           data = pp_nomag_dat(),
                           init = list(list(theta = 0.0)),
                           mcmc_args = list(chains = 2L, iter_warmup = 100L,
                                             iter_sampling = 100L, seed = 7L,
                                             refresh = 0, show_messages = FALSE)),
    "magnitude_adj")
})

test_that("refit_magnitude_adjust validates the refit model, not just the fit", {
  skip_stan()
  # The fit HAS magnitude_adj (truelik) but the model being refit LACKS it
  # (nomag). A fit-only check would pass here; the model-side check must catch
  # the mismatch, since the rescale is applied through `model`.
  fit <- get_truelik_fit()
  expect_error(
    refit_magnitude_adjust(fit, model = get_nomag_model(),
                           data = pp_nomag_dat(),
                           init = list(list(theta = 0.0)),
                           mcmc_args = list(chains = 2L, iter_warmup = 100L,
                                             iter_sampling = 100L, seed = 7L,
                                             refresh = 0, show_messages = FALSE)),
    "magnitude_adj")
})

test_that("check_model_data_arg flags a model lacking the data argument", {
  skip_stan()
  # nomag lacks magnitude_adj; truelik declares it.
  expect_error(pseudopost:::check_model_data_arg(get_nomag_model(), "magnitude_adj"),
               "magnitude_adj")
  expect_no_error(pseudopost:::check_model_data_arg(get_truelik_model(), "magnitude_adj"))
})

test_that("check_pl_requirements(need_magnitude = TRUE) flags the nomag fit but not truelik", {
  skip_stan()
  # nomag has log_lik but no magnitude_adj -> the magnitude requirement fails.
  expect_error(check_pl_requirements(get_nomag_fit(), need_magnitude = TRUE),
               "magnitude_adj")
  # truelik has magnitude_adj -> passes.
  expect_no_error(check_pl_requirements(get_truelik_fit(), need_magnitude = TRUE))
})
