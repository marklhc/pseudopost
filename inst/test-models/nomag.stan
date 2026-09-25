//
// Minimal per-observation likelihood fixture that deliberately OMITS the
// `magnitude_adj` data argument (and `use_priors`). Used to verify that the
// magnitude-adjustment entry points reject a model that cannot be rescaled.
//
data {
  int<lower=1> N;
  vector[N] y;
}
parameters {
  real theta;
}
model {
  for (n in 1:N) target += normal_lpdf(y[n] | theta, 1);
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N) log_lik[n] = normal_lpdf(y[n] | theta, 1);
}
