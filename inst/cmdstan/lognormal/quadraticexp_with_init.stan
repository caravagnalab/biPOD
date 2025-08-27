functions {
  real integrated_r(real t, real t0, vector t_array, vector rho_array) {
    int n_t = num_elements(t_array);
    real res = 0;

    if (n_t == 0) return rho_array[1] * (t - t0);
    if (t <= t_array[1]) return rho_array[1] * (t - t0);

    res = rho_array[1] * (t_array[1] - t0);
    for (i in 2:n_t) {
      if (t <= t_array[i])
        return res + rho_array[i] * (t - t_array[i - 1]);
      res += rho_array[i] * (t_array[i] - t_array[i - 1]);
    }
    return res + rho_array[n_t + 1] * (t - t_array[n_t]);
  }

  real mean_t(real t, real t0, vector t_array, vector rho_array, real b) {
    real dt = t - t0;
    return exp( integrated_r(t, t0, t_array, rho_array) + b * dt * dt );
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
  vector[G] rho;
  real<upper=T[1]> t0;
  real b;                    // quadratic curvature
  real<lower=0> sigma;
}

model {
  // Priors (adjust as needed)
  rho ~ normal(0, 1);
  t0  ~ normal(T[1], 30);
  b   ~ normal(0, 0.5);      // modest curvature prior
  sigma ~ exponential(1);

  if (prior_only == 0) {
    vector[S] mu_log; // location for LogNormal (mean on log scale)
    for (i in 1:S) {
      real mu_pred = mean_t(T[i], t0, t_array, rho, b);
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
      mu_pred[i] = mean_t(T[i], t0, t_array, rho, b);
      log_lik[i] = lognormal_lpdf(N[i] | log(mu_pred[i]), sigma);
      yrep[i]    = lognormal_rng(log(mu_pred[i]), sigma);
    }
  }
}
