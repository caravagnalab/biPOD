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

  for (i in 1:S) {
    real mu = mean_t(T[i], t0, t_array, rho, b);
    mu = fmax(mu, 1e-6);
    log(N[i] + 1e-3) ~ normal(log(mu), sigma);
  }
}

generated quantities {
  vector[S] log_lik;
  vector[S] yrep;
  vector[S] mu_pred;

  for (i in 1:S) {
    mu_pred[i] = mean_t(T[i], t0, t_array, rho, b);
    mu_pred[i] = fmax(mu_pred[i], 1e-6);
    yrep[i] = exp(normal_rng(log(mu_pred[i]), sigma));
    log_lik[i] = normal_lpdf(log(N[i] + 1e-3) | log(mu_pred[i]), sigma);
  }
}
