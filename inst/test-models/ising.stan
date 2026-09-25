//
// Besag pseudo-likelihood for a fully-connected binary Ising model
// (e.g., a PTSD symptom network). Common coupling beta;
// symptom-specific intercepts alpha[i]. log_lik (one entry per
// individual) computed in generated quantities to avoid autodiff
// overhead during sampling.
//
data {
  int<lower=1> n;
  int<lower=2> p;
  array[n, p] int<lower=0, upper=1> x;
  real<lower=0> magnitude_adj;
  int<lower=0, upper=1> use_priors;
}
transformed data {
  array[p, n] int<lower=0, upper=1> xt;
  matrix[n, p] S;
  {
    for (r in 1:n) {
      real total_active = sum(x[r, ]);
      for (i in 1:p) {
        xt[i, r] = x[r, i];
        S[r, i] = total_active - x[r, i];
      }
    }
  }
}
parameters {
  vector[p] alpha;
  real beta;
}
model {
  if (use_priors == 1) {
    alpha ~ std_normal();
    beta ~ std_normal();
  }
  {
    real pll = 0;
    for (i in 1:p) {
      pll += bernoulli_logit_lpmf(xt[i] | alpha[i] + beta * S[, i]);
    }
    target += magnitude_adj * pll;
  }
}
generated quantities {
  vector[n] log_lik = rep_vector(0.0, n);
  {
    for (i in 1:p) {
      vector[n] eta_i = alpha[i] + beta * S[, i];
      for (r in 1:n) {
        log_lik[r] += bernoulli_logit_lpmf(xt[i, r] | eta_i[r]);
      }
    }
  }
}
