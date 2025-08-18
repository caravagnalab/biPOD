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

    // Add numerical stability checks
    if (n0 <= 0 || K <= 0) return 1e-6;
    if (K <= n0) return fmin(n0, K);

    real ratio = (K - n0) / n0;
    real exp_arg = -rint;

    if (exp_arg > 700) return n0;
    if (exp_arg < -700) return K;

    real exp_term = exp(exp_arg);
    real denominator = 1 + ratio * exp_term;

    if (denominator <= 0) return 1e-6;

    real result = K / denominator;
    return (is_nan(result) || result <= 0) ? 1e-6 : result;
  }
}

data {
  int<lower=1> S;
  int<lower=1> G;
  array[S] int<lower=0> N;
  array[S] real T;
  vector[G - 1] t_array;
}

parameters {
  real<lower=0> n0;
  real<lower=0> K;
  vector[G] rho;
  real<lower=0> sigma;  // Log-normal noise
}

transformed parameters {
  real K_constrained = n0 + K;
}

model {
  rho ~ normal(0, 1);
  n0 ~ normal(N[1], N[1] / 20.0);
  K ~ lognormal(log(max(N) + 1), 1.0);
  sigma ~ exponential(1);

  for (i in 1:S) {
    real mu = mean_t(T[i], T[1], n0, K_constrained, t_array, rho);
    mu = fmax(mu, 1e-6);
    log(N[i] + 1e-3) ~ normal(log(mu), sigma);
  }
}

generated quantities {
  vector[S] log_lik;
  vector[S] yrep;
  vector[S] mu_pred;

  for (i in 1:S) {
    mu_pred[i] = mean_t(T[i], T[1], n0, K_constrained, t_array, rho);
    mu_pred[i] = fmax(mu_pred[i], 1e-6);
    yrep[i] = exp(normal_rng(log(mu_pred[i]), sigma));
    log_lik[i] = normal_lpdf(log(N[i] + 1e-3) | log(mu_pred[i]), sigma);
  }
}
