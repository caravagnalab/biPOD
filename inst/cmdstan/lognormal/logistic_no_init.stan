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

  // Logistic mean with piecewise r(t); numerically stable
  real mean_t(real t, real t0, real n0, real K, vector t_array, vector rho_array) {
    // Preconditions guard
    if (n0 <= 0 || K <= 0) return 1e-6;

    real rint = integrated_r(t, t0, t_array, rho_array);
    // y(t) = K / (1 + ((K-n0)/n0) * exp(-rint))
    // Use exp over/underflow guards
    real log_ratio = log_diff_exp(log(K), log(n0)) - log(n0); // log((K-n0)/n0)
    real z = log1p_exp(log_ratio - rint); // log(1 + ((K-n0)/n0)*exp(-rint))
    // y = K / (1 + ...) = exp(log(K) - z)
    real y = exp(log(K) - z);
    return fmax(y, 1e-6);
  }
}

data {
  int<lower=1> S;
  int<lower=1> G;
  array[S] int<lower=1> N;
  array[S] real T;
  vector[G - 1] t_array;
  int<lower=0,upper=1> prior_only;
}

parameters {
  real<lower=1> K;
  vector[G] rho;
  real<lower=0> n0;
  real<lower=0> sigma;  // for log-normal observation noise
}

model {
  rho ~ normal(0, 1);
  n0 ~ normal(N[1], N[1] / 5.0);
  K ~ normal(max(N), max(N));
  sigma ~ exponential(1);

  if (prior_only == 0) {
    vector[S] mu_log; // location for LogNormal (mean on log scale)
    for (i in 1:S) {
      real mu_pred = mean_t(T[i], T[1], n0, K, t_array, rho);
      mu_log[i] = log(mu_pred);
    }
    N ~ lognormal(mu_log, sigma);
  }
}

generated quantities {
  vector[S] log_lik;
  array[S] real yrep;
  vector[S] mu_pred;

  if (prior_only == 0) {
    for (i in 1:S) {
      mu_pred[i] = mean_t(T[i], T[1], n0, K, t_array, rho);
      log_lik[i] = lognormal_lpdf(N[i] | log(mu_pred[i]), sigma);
      yrep[i]    = lognormal_rng(log(mu_pred[i]), sigma);
    }
  }
}
