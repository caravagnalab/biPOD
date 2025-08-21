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
    // y(t) = K - (K - n0) * exp(-∫ rho)
    return K - (K - n0) * exp(-integrated_r(t, t0, t_array, rho_array));
  }
}

data {
  int<lower=1> S;
  int<lower=1> G;                 // length(rho) = G
  array[S] int<lower=0> N;
  array[S] real T;
  vector[G - 1] t_array;          // breakpoints (increasing)
}

parameters {
  real<lower=0> n0;               // intercept at t = T[1]
  vector[G] rho;                  // segment rates
  real<lower=0> K_plus;           // ensures K > n0 via K = n0 + exp(K_plus)
  real<lower=0> sigma;            // log-normal noise sd
}

transformed parameters {
  real K = n0 + exp(K_plus);      // enforce K > n0 without parameter-dependent bounds
}

model {
  // Priors (tune as needed)
  rho ~ normal(0, 1);
  n0 ~ normal(N[1], N[1] / 20.0);
  K_plus ~ normal(0, 1);
  sigma ~ exponential(1);

  for (i in 2:S) {
    real mu = mean_t(T[i], T[1], n0, K, t_array, rho);
    mu = fmax(mu, 1e-6);
    log(N[i] + 1e-3) ~ normal(log(mu), sigma);
  }
}

generated quantities {
  vector[S] log_lik;
  vector[S] yrep;
  vector[S] mu_pred;

  for (i in 1:S) {
    mu_pred[i] = mean_t(T[i], T[1], n0, K, t_array, rho);
    mu_pred[i] = fmax(mu_pred[i], 1e-6);
    yrep[i] = exp(normal_rng(log(mu_pred[i]), sigma));
    log_lik[i] = (i >= 2)
      ? normal_lpdf(log(N[i] + 1e-3) | log(mu_pred[i]), sigma)
      : 0;
  }
}
