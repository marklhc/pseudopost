//
// Estimates a 1-factor CFA model based on pairwise polychoric composite likelihood.
// (Theta Parameterization for HMC Efficiency)
//
#include "pairwise_fun.stan"
data {
  int<lower=1> p;
  int<lower=1> k;

  int<lower=1> num_pairs;
  array[num_pairs, k, k] int<lower=0> ns;

  int<lower=1> N;
  array[N, p] int<lower=1, upper=k> y;

  real<lower=0> magnitude_adj;
  int<lower=0, upper=1> use_priors;
}
transformed data {
  array[num_pairs, 2] int pair_idx;
  {
    int s = 1;
    for (i in 1:(p - 1)) {
      for (j in (i + 1):p) {
        pair_idx[s, 1] = i;
        pair_idx[s, 2] = j;
        s += 1;
      }
    }
  }

  array[p, N] int yt;
  for (n in 1:N)
    for (i in 1:p)
      yt[i, n] = y[n, i];
}
parameters {
  real<lower=0> lambda_1;
  vector[p - 1] lambda_rest;
  array[p] ordered[k - 1] thres;
}

transformed parameters {
  vector[p] lambda = append_row(lambda_1, lambda_rest);
  vector[p] lambda_delta;
  array[p] vector[k - 1] thres_delta;
  array[p] vector[k - 1] phi_thres;

  for (i in 1:p) {
    real sd_total = sqrt(square(lambda[i]) + 1.0);
    lambda_delta[i] = lambda[i] / sd_total;
    thres_delta[i] = thres[i] / sd_total;
    phi_thres[i] = Phi(thres_delta[i]);
  }
}

model {
  if (use_priors == 1) {
    lambda ~ normal(0.5, 1);
    for (i in 1:p) {
      thres[i] ~ normal(0, 3);
    }
  }

  {
    real pll = 0.0;
    for (s in 1:num_pairs) {
      int pi = pair_idx[s, 1];
      int pj = pair_idx[s, 2];
      matrix[k, k] logprob =
        polycor_to_logprob(thres_delta[pi], thres_delta[pj],
                           phi_thres[pi], phi_thres[pj],
                           lambda_delta[pi] * lambda_delta[pj]);

      for (c1 in 1:k) {
        for (c2 in 1:k) {
          if (ns[s, c1, c2] > 0) {
            pll += ns[s, c1, c2] * logprob[c1, c2];
          }
        }
      }
    }
    target += magnitude_adj * pll;
  }
}

generated quantities {
  vector[N] log_lik = rep_vector(0.0, N);

  {
    for (s in 1:num_pairs) {
      int pi = pair_idx[s, 1];
      int pj = pair_idx[s, 2];
      matrix[k, k] logprob =
        polycor_to_logprob(thres_delta[pi], thres_delta[pj],
                           phi_thres[pi], phi_thres[pj],
                           lambda_delta[pi] * lambda_delta[pj]);

      for (n in 1:N) {
        log_lik[n] += logprob[yt[pi, n], yt[pj, n]];
      }
    }
  }
}
