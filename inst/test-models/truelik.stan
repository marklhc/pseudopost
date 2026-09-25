data {
  int<lower=1> N;
  vector[N] y;
  real<lower=0> magnitude_adj;
  int<lower=0, upper=1> use_priors;
}
parameters { real theta; }
model {
  if (use_priors == 1) theta ~ normal(0, 10);
  for (n in 1:N) target += magnitude_adj * normal_lpdf(y[n] | theta, 1);
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N) log_lik[n] = normal_lpdf(y[n] | theta, 1);
}
