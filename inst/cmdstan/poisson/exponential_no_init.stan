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

  real mean_t(real t, real t0, real n0, vector t_array, vector rho_array) {
    return n0 * exp(integrated_r(t, t0, t_array, rho_array));
  }
}

data {
  int<lower=1> S;
  int<lower=1> G;
  array[S] int<lower=0> N;
  array[S] real T;
  vector[G - 1] t_array;
  int<lower=0,upper=1> prior_only;
}

parameters {
  real<lower=0> n0;
  vector[G] rho;
}

model {
  rho ~ normal(0, 1);
  n0 ~ normal(N[1], N[1] / 5.0);
  

  if (prior_only == 0) {
    vector[S] mu_pred;
    for (i in 1:S) {
      mu_pred[i] = mean_t(T[i], T[1], n0, t_array, rho);
    }
    N ~ poisson(mu_pred);
  }
}


generated quantities {
  vector[S] log_lik;
  array[S] real yrep;
  vector[S] mu_pred;

  if (prior_only == 0) {
    for (i in 1:S) {
      mu_pred[i] = mean_t(T[i], T[1], n0, t_array, rho);
      log_lik[i] = poisson_lpmf(N[i] | mu_pred[i]);
      yrep[i]    = poisson_rng(mu_pred[i]);
    }
  }
}

