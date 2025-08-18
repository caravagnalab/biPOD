functions {
  real integrated_r(real t, real t0, vector t_array, vector rho_array) {
    int n_t = num_elements(t_array);
    int n_rho = num_elements(rho_array);
    real res = 0;
    if (n_t == 0) return rho_array[1] * (t - t0);
    if (t <= t_array[1]) return rho_array[1] * (t - t0);
    res = rho_array[1] * (t_array[1] - t0);
    for (i in 2:n_t) {
      if (t <= t_array[i])
        return res + rho_array[i] * (t - t_array[i - 1]);
      res += rho_array[i] * (t_array[i] - t_array[i - 1]);
    }
    res += rho_array[n_rho] * (t - t_array[n_t]);
    return res;
  }

  real mean_t(real t, real t0, real n0, real K, vector t_array, vector rho_array) {
    real rint = integrated_r(t, t0, t_array, rho_array);

    // Very conservative numerical protection
    if (K <= 1e-6 || n0 <= 1e-10) return 1e-6;

    // Ensure n0/K is well within (0,1)
    real ratio = n0 / K;
    ratio = fmax(1e-8, fmin(1 - 1e-8, ratio));

    // Compute C = log(-log(ratio)) with protection
    real log_ratio = log(ratio);
    real neg_log_ratio = -log_ratio;
    real C = log(neg_log_ratio);

    // Compute the final result with clipping
    real exp_inner = fmax(-50, fmin(50, -rint + C));  // Tighter bounds
    real exp_outer = fmax(-50, fmin(50, -exp(exp_inner)));

    real result = K * exp(exp_outer);
    return fmax(1e-6, fmin(result, K + 1e-6));
  }
}

data {
  int<lower=1> S;
  int<lower=1> G;
  array[S] int<lower=0> N;
  array[S] real T;
  vector[G - 1] t_array;
}

transformed data {
  // Compute some useful quantities for priors
  real min_N = min(N);
  real max_N = max(N);
  real mean_N = mean(to_vector(N));
  real sd_N = sd(to_vector(N));
}

parameters {
  real<lower=0> K_raw;        // Raw parameter for K
  real<lower=0> n0_ratio;     // n0 as fraction of K
  vector[G] rho_raw;          // Raw parameters for rho
  real<lower=0> sigma;  // standard deviation of log-observation noise
}

transformed parameters {
  // More stable parameterization
  real K = mean_N + K_raw * (max_N - mean_N + sd_N);  // K around observed range
  real n0 = n0_ratio * K * 0.95;  // n0 is fraction of K, max 95% of K
  vector[G] rho = rho_raw * 0.1;  // Scale down rho to prevent extreme values
}

model {
  // Informative priors based on data
  K_raw ~ exponential(2);     // Keeps K reasonable
  n0_ratio ~ beta(2, 3);      // Favors n0 being smaller fraction of K
  rho_raw ~ std_normal();     // Standard normal, scaled in transformed params
  sigma ~ exponential(1);

  for (i in 1:S) {
    real mu = mean_t(T[i], T[1], n0, K, t_array, rho);
    mu = fmax(mu, 1e-6);  // still guard against log(0)
    log(N[i] + 1e-3) ~ normal(log(mu), sigma);  // log-normal likelihood
  }
}

generated quantities {
  vector[S] log_lik;
  vector[S] yrep;
  vector[S] mu_pred;

  for (i in 1:S) {
    mu_pred[i] = mean_t(T[i], T[1], n0, K, t_array, rho);
    mu_pred[i] = fmax(mu_pred[i], 1e-6);
    
    // Sample from the log-normal predictive distribution
    yrep[i] = exp(normal_rng(log(mu_pred[i]), sigma));
    
    // Compute log-likelihood using the same log-normal model
    log_lik[i] = normal_lpdf(log(N[i] + 1e-3) | log(mu_pred[i]), sigma);
  }
}
