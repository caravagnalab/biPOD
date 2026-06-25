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
    real log_ratio = log(abs(K - n0)) - log(n0);  // fabs -> abs
    real z = log_ratio - rint;
    return exp(log(K) - log1p_exp(z));
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
transformed data {
  real max_N = max(N);
}
parameters {
  real<lower=max_N * 0.9> K;  // carrying capacity, at least 90% of max observed
  vector[G] rho;
  real<lower=0> n0;
  real<lower=0> phi;
}
model {
  rho ~ normal(0, 1);
  n0  ~ normal(N[1], N[1] / 5.0);
  K   ~ normal(max_N, max_N);
  phi ~ gamma(100, 1);
  if (prior_only == 0) {
    vector[S] mu_pred;
    for (i in 1:S) {
      mu_pred[i] = mean_t(T[i], T[1], n0, K, t_array, rho);
    }
    N ~ neg_binomial_2(mu_pred, phi);
  }
}
generated quantities {
  vector[S] log_lik;
  array[S] int yrep;
  vector[S] mu_pred;
  if (prior_only == 0) {
    for (i in 1:S) {
      mu_pred[i] = mean_t(T[i], T[1], n0, K, t_array, rho);
      log_lik[i] = neg_binomial_2_lpmf(N[i] | mu_pred[i], phi);
      yrep[i]    = neg_binomial_2_rng(mu_pred[i], phi);
    }
  }
}
