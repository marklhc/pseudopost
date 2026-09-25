functions {
  
  real binormal_cdf(tuple(real, real) z, real rho, tuple(real, real) phi_z) {
    // Boundary protection for mathematical limits
    if (rho <= -1.0 || rho >= 1.0) reject("binormal_cdf bounds error: rho = ", rho);
    if (rho == 1.0) return min([phi_z.1, phi_z.2]); 
    if (rho == -1.0) return fmax(0.0, phi_z.1 - (1.0 - phi_z.2)); 
    if (z.1 == 0 && z.2 == 0) return 0.25 + asin(rho) / (2 * pi());
    
    real denom = sqrt((1 + rho) * (1 - rho));
    real term1 = z.1 == 0
      ? (z.2 > 0 ? 0.25 : -0.25)
      : owens_t(z.1, (z.2 / z.1 - rho) / denom);
    real term2 = z.2 == 0
      ? (z.1 > 0 ? 0.25 : -0.25)
      : owens_t(z.2, (z.1 / z.2 - rho) / denom);
    real z1z2 = z.1 * z.2;
    real delta = z1z2 < 0 || (z1z2 == 0 && (z.1 + z.2) < 0);
    return 0.5 * (phi_z.1 + phi_z.2 - delta) - term1 - term2;
  }
  
  // Returns log-probabilities directly, utilizing the zero-padding boundary trick
  matrix polycor_to_logprob(vector b1, vector b2, vector phi_b1, vector phi_b2, real rho) {
    int k1 = size(b1) + 1;
    int k2 = size(b2) + 1;
    
    matrix[k1 + 1, k2 + 1] cprob = rep_matrix(0.0, k1 + 1, k2 + 1);
    matrix[k1, k2] raw_prob;
    
    // Fill margins and corner
    for (i in 1:(k1 - 1)) cprob[i + 1, k2 + 1] = phi_b1[i];
    for (j in 1:(k2 - 1)) cprob[k1 + 1, j + 1] = phi_b2[j];
    cprob[k1 + 1, k2 + 1] = 1.0;
    
    // Fill interior
    for (i in 1:(k1 - 1)) {
      for (j in 1:(k2 - 1)) {
        cprob[i + 1, j + 1] = binormal_cdf((b1[i], b2[j]) | rho, (phi_b1[i], phi_b2[j]));
      }
    }
    
    // Unified probability calculation with underflow protection
    real total_prob = 0.0;
    for (i in 1:k1) {
      for (j in 1:k2) {
        real p = cprob[i + 1, j + 1] - cprob[i, j + 1] - cprob[i + 1, j] + cprob[i, j];
        raw_prob[i, j] = fmax(1e-10, p);
        total_prob += raw_prob[i, j];
      }
    }
    
    return log(raw_prob / total_prob);
  }
  
  // Computes the composite log-likelihood without the multinomial constant
  real pll_1fcfa(array[,,] int ns, array[] vector thres, array[] vector phi_thres, vector lambda) {
    int p = size(lambda);
    real out = 0.0;
    int s = 1;
    
    for (i in 1:(p - 1)) {
      for (j in (i + 1):p) {
        matrix[size(thres[i]) + 1, size(thres[j]) + 1] logprob = 
          polycor_to_logprob(thres[i], thres[j], phi_thres[i], phi_thres[j], lambda[i] * lambda[j]);
          
        for (c1 in 1:(size(thres[i]) + 1)) {
          for (c2 in 1:(size(thres[j]) + 1)) {
            // Bypass empty cells entirely
            if (ns[s, c1, c2] > 0) {
              out += ns[s, c1, c2] * logprob[c1, c2];
            }
          }
        }
        s += 1;
      }
    }
    return out;
  }

  // Optimized for LOO using direct matrix lookups
  vector pll_1fcfa_is(array[,] int y, array[] vector thres, vector lambda) {
    int N = size(y);
    int p = size(lambda);
    int k = size(thres[1]) + 1; 
    
    vector[N] out = rep_vector(0.0, N);
    array[p] vector[k - 1] phi_thres;
    
    for (i in 1:p) {
      phi_thres[i] = Phi(thres[i]);
    }
    
    for (i in 1:(p - 1)) {
      for (j in (i + 1):p) {
        // Fetch the fully calculated log-probability grid
        matrix[k, k] logprob = 
          polycor_to_logprob(thres[i], thres[j], phi_thres[i], phi_thres[j], lambda[i] * lambda[j]);
          
        for (n in 1:N) {
          // Instant lookup; zero mathematical operations inside the N loop
          out[n] += logprob[y[n, i], y[n, j]]; 
        }
      }
    }
    return out;
  }
}