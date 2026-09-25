//
// Exact likelihood for a fully-connected Ising model with common coupling.
// The partition function is computed exactly by summing over all 2^p configs.
//
data {
  int<lower=1> n;
  int<lower=2> p;
  array[n, p] int<lower=0, upper=1> x;
  int<lower=1> num_configs;
  matrix[num_configs, p] all_configs;
}
parameters {
  vector[p] alpha;
  real beta;
}
transformed parameters {
  vector[num_configs] config_energies;
  real logZ;
  {
    for (c in 1:num_configs) {
      real e = 0;
      for (i in 1:p) e += alpha[i] * all_configs[c, i];
      for (i in 1:(p - 1))
        for (j in (i + 1):p)
          e += beta * all_configs[c, i] * all_configs[c, j];
      config_energies[c] = e;
    }
    logZ = log_sum_exp(config_energies);
  }
}
model {
  alpha ~ std_normal();
  beta ~ std_normal();
  for (r in 1:n) {
    real e = 0;
    for (i in 1:p) e += alpha[i] * x[r, i];
    for (i in 1:(p - 1))
      for (j in (i + 1):p)
        e += beta * x[r, i] * x[r, j];
    target += e - logZ;
  }
}
