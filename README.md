
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pseudopost

**Post-hoc pseudo-posterior uncertainty adjustment for composite and
pseudo-likelihood Stan models.**

When a Stan model is fit to a *composite*, *pairwise*, or otherwise
*pseudo-likelihood* target (via
[cmdstanr](https://mc-stan.org/cmdstanr/)), the posterior that MCMC
recovers is the posterior of the pseudo-likelihood. For such models the
**naive** posterior covariance — the inverse Hessian, or simply the
sample covariance of the draws — is generally *miscalibrated*: standard
errors and credible intervals are often, but not always, too small and
too narrow (i.e. *anti-conservative*). `pseudopost` provides post-hoc
covariance (infinitesimal jackknife / Godambe sandwich) adjustments that
target the robust covariance $V = H^{-1} J H^{-1}$ to improve scale and
coverage **post-hoc**, without refitting the model — an approximation to
the correct uncertainty, not a guarantee of correct coverage.

## What it does

`pseudopost` implements post-hoc covariance (sandwich) adjustments, all
driven by a single entry point, `adjust_pseudo_posterior()`, organized
around the robust covariance $V = H^{-1} J H^{-1}$ (Shaby, 2014):

- an **infinitesimal jackknife** (`"ij"`) estimator — a model-blind
  estimate of the full $V$ straight from the MCMC draws (Giordano &
  Broderick, 2023);
- a **closed** (`"sandwich"`) Godambe estimator with an explicit
  precision `bread` ($H$) and score `meat` ($J$);
- an **open-faced** (`"ofs"`) sandwich estimator (Shaby, 2014): the same
  explicit `bread` ($Q$) and `meat` ($P$), applied as the draw transform
  $\Omega = (Q^{-1} P)^{1/2}$ without calibrating to the empirical draw
  covariance (the “one-bread” assumption);
- a `ci` type (percentile of the adjusted draws, or a Wald interval).

The four axes are orthogonal and each may be a vector, so they run a
sub-grid (one table row per parameter and valid combination).

| axis | values | role |
|----|----|----|
| `estimator` | `"ij"`, `"sandwich"`, `"ofs"` | infinitesimal jackknife (model-blind $V$ from the draws), closed Godambe $H^{-1} J H^{-1}$, or open-faced draw transform $(Q^{-1} P)^{1/2}$ |
| `bread` | `"mcmc"`, `"hessian_mle"`, `"hessian_map"` | the precision $H$; `"sandwich"`/`"ofs"` only |
| `meat` | `"score"`, `"score_avg"` | the score moment $J$; `"sandwich"`/`"ofs"` only |
| `ci` | `"quantile"`, `"wald"` | percentile of the adjusted draws vs estimate $\pm z \cdot \text{se}$ |

The `"ij"` estimator does not use `bread`/`meat` (only `ci` varies). The
`"sandwich"` and `"ofs"` combinations are the `bread` × `meat` pairs
below, each also crossed with `ci`:

| `bread`         | `meat`        |
|-----------------|---------------|
| `"mcmc"`        | `"score"`     |
| `"mcmc"`        | `"score_avg"` |
| `"hessian_mle"` | `"score"`     |
| `"hessian_mle"` | `"score_avg"` |
| `"hessian_map"` | `"score"`     |
| `"hessian_map"` | `"score_avg"` |

The baseline (no robust covariance) is the legacy `methods = "noadj"`
preset.

A scalar **magnitude** (scale) adjustment — the Godambe ratio
$m = d / \text{tr}(H^{-1} J)$ — is available via `magnitude_adjust()`
(compute $m$) and `refit_magnitude_adjust()` (compute $m$, then refit
with `magnitude_adj = m`). This is a single *global* scalar that
rescales the whole pseudo-likelihood: a coarse correction of the overall
posterior scale that, in general, does not reproduce a
direction-dependent multivariate sandwich covariance.

> **See** *The pseudo-likelihood requirements* (vignette,
> `vignettes/pseudo-likelihood-requirements.Rmd`) for what a qualifying
> model must expose and a full worked walkthrough.

## Installation

``` r
# From GitHub
remotes::install_github("marklai/pseudopost")
```

`pseudopost` depends on [`cmdstanr`](https://mc-stan.org/cmdstanr/) and
a working [CmdStan](https://mc-stan.org/cmdstan/) toolchain to inspect
and (re)optimise fitted models.

## Quick start

The package ships a one-parameter normal model (`truelik.stan`) that is
a *true* per-observation likelihood, useful as a sanity check: every
correction must collapse to the identity.

``` r
library(cmdstanr)
library(pseudopost)

model <- cmdstan_model(
  system.file("test-models", "truelik.stan", package = "pseudopost"),
  force_recompile = TRUE, quiet = TRUE)
set.seed(42)
dat <- list(N = 100L, y = rnorm(100, 1, 1),
            magnitude_adj = 1.0, use_priors = 1L)
fit <- model$sample(data = dat, chains = 2L, iter_warmup = 200L,
                    iter_sampling = 500L, seed = 42L, refresh = 0,
                    show_messages = FALSE)

# Validate the model meets the requirements
check_pl_requirements(fit)

# Run the adjustment grid
pp <- adjust_pseudo_posterior(
  fit, params = "theta",
  estimator = c("ij", "sandwich"),
  bread = c("mcmc", "hessian_map"),
  meat = "score",
  ci = "quantile",
  model = model, data = dat, init = list(list(theta = 0.0)),
  verbose = FALSE)
print(pp)
```

For your own model, the workflow is the same:

``` r
model <- cmdstan_model("my_pseudo_model.stan")
fit   <- model$sample(data = my_dat)

# The model must expose a per-observation log-likelihood generated quantity
# `log_lik` (one entry per pseudo-observation); see the requirements vignette.
check_pl_requirements(fit)

pp <- adjust_pseudo_posterior(
  fit, params = c("lambda", "thres"),
  estimator = "sandwich",
  bread = c("mcmc", "hessian_map"),
  meat = "score",
  ci = "quantile",
  model = model, data = my_dat)
print(pp)
```

`pp$table` is a long-form data frame (`method, param, est, se, lo, hi`);
`pp$draws` holds the affine-adjusted posterior draws per method;
`pp$meta` carries the fit, settings, and any warnings/notes.

## The adjustment families

All covariance adjustments estimate the robust covariance
$V = H^{-1} J H^{-1}$ (Shaby, 2014) and differ in *how* that
decomposition is estimated:

- **Infinitesimal jackknife (`estimator = "ij"`)** — a model-blind
  estimate of the full $V$ straight from the MCMC draws (Giordano &
  Broderick, 2023); the CRAN
  [`IJSE`](https://cran.r-project.org/package=IJSE) package (Ji, Lee &
  Rabe-Hesketh, 2025) implements the same idea. It does not use
  `bread`/`meat`.
- **Closed Godambe sandwich (`estimator = "sandwich"`)** — the explicit
  $V = H^{-1} J H^{-1}$, with:
  - `bread` ($H$, precision): `"mcmc"` (inverse MCMC draw covariance),
    `"hessian_mle"` (negative Hessian of the pseudo-likelihood at the
    MLE, no Jacobian), or `"hessian_map"` (negative Hessian of the
    pseudo-posterior at the MAP, with Jacobian);
  - `meat` ($J$): `"score"` (the moment estimator from the
    per-observation score vector) or `"score_avg"` (a draw-averaged
    version).
- **Open-faced sandwich (`estimator = "ofs"`)** — Shaby’s (2014)
  original adjustment: the same explicit `bread` ($Q$) and `meat` ($P$)
  as the closed sandwich, but applied directly to the centered draws as
  $\Omega = (Q^{-1} P)^{1/2}$ without calibrating to the empirical draw
  covariance $C_{\text{emp}}$ (the “one-bread” assumption
  $C_{\text{emp}} = Q^{-1}$).
- **Magnitude** — a single *global* scalar that rescales the whole
  pseudo-likelihood; a coarse correction of the overall posterior scale
  that, in general, does not reproduce a direction-dependent
  multivariate sandwich covariance (see above). \<!– The original named
  presets remain as **silent aliases** for the axes above, so existing
  code keeps working; when `methods` is supplied it takes precedence
  over `estimator`/`bread`/`meat`/`ci`:

| legacy `methods` | `estimator`  | `bread`         | `meat`    | `ci`         |
|------------------|--------------|-----------------|-----------|--------------|
| `"noadj"`        | (baseline)   | —               | —         | —            |
| `"ijse"`         | `"ij"`       | —               | —         | `"wald"`     |
| `"ijse2"`        | `"ij"`       | —               | —         | `"quantile"` |
| `"curvadj"`      | `"sandwich"` | `"hessian_mle"` | `"score"` | `"quantile"` |
| `"curvadj2"`     | `"sandwich"` | `"hessian_map"` | `"score"` | `"quantile"` |
| `"curvadj3"`     | `"sandwich"` | `"mcmc"`        | `"score"` | `"quantile"` |

## Scope and limits

- The `"hessian_mle"` / `"hessian_map"` breads evaluate the Hessian at
  an MLE/MAP, so they require a `model`/`data` pair (or a pre-optimized
  `opt`/`opt_map`); there is **no MCMC-posterior-mean fallback** (the
  posterior mean is not an MLE or MAP and would mislabel the point). For
  `"hessian_mle"` to be a true likelihood-only MLE, the model must
  disable **all** priors when `use_priors = 0` — the package cannot
  verify that the Stan program honours that switch.
- The affine draw transformation sets the adjusted draws’ covariance to
  the target in the **unconstrained** space (the map is built from the
  empirical draw covariance); mapping back to the constrained `params`
  can change the covariance in the constrained space, since the affine
  map is not preserved under non-linear transforms.
- The `log_lik` requirement is a **vector** of per-observation
  contributions (one entry per pseudo-observation); the package checks
  the shape only and does not verify that the entries are contributions
  to independent observational units or that they sum to the fitted
  pseudo-likelihood.

## References

- Shaby, B. A. (2014). The Open-Faced Sandwich Adjustment for MCMC Using
  Estimating Functions. *Journal of Computational and Graphical
  Statistics*, 23(3), 853-876. <doi:10.1080/10618600.2013.842174>.
- Giordano, R., & Broderick, T. (2023). The Bayesian Infinitesimal
  Jackknife for Variance. arXiv preprint arXiv:2305.06466.
  <doi:10.48550/arXiv.2305.06466>.
- Ji, F., Lee, J., & Rabe-Hesketh, S. (2025). Valid standard errors for
  Bayesian quantile regression with clustered and independent data.
  *Journal of Educational and Behavioral Statistics*, 51(5), 1055-1084.
  <doi:10.3102/10769986251379738>.
- Godambe, V. P., & Heyde, C. C. (1975). Quasi-likelihood and optimal
  estimation. *Journal of the American Statistical Association*,
  70(351), 455-467.

## License

GPL (\>= 3)
