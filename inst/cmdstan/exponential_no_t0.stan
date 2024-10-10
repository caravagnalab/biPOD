
functions {
  real mean_t (real t, real t0, real n0, vector t_array, vector rho_array) {
    real res;
    int n_t = num_elements(t_array);
    int n_rho = num_elements(rho_array);

    if (n_t == 0) return n0 * exp(rho_array[1] * (t - t0));

    if (t <= t_array[1]) return n0 * exp(rho_array[1] * (t - t0));

    res = n0 * exp(rho_array[1] * (t_array[1] - t0));

    for (i in 2:n_t) {
      if (t <= t_array[i]) {
        return res * exp(rho_array[i] * (t - t_array[i-1]));
      } else {
        res = res * exp(rho_array[i] * (t_array[i] - t_array[i-1]));
      }
    }
    res = res * exp(rho_array[n_rho] * (t - t_array[n_t]));
    return res;
  }
}

data {
  int<lower=1> S; // Number of steps
  int<lower=1> G; // Number of windows

  array[S] int<lower=0> N; // observations
  array[S] real T;         // observations

  vector[G-1] t_array;
}

parameters {
  real n0;
  vector[G] rho;
}

model {
  for (i in 1:G) {
    target += normal_lpdf(rho[i] | 0, 1);
  }

  target += normal_lpdf(n0 | N[1], N[1] / 20.0);

  for (i in 2:S) {
    target += poisson_lpmf(N[i] | mean_t(T[i], T[1], n0, t_array, rho));
  }
}

generated quantities {
  vector[S] log_lik;             // Log-likelihood for each observation
  vector[S] yrep;               // Expected values for N given x

  for (i in 1:S) {
    yrep[i] = mean_t(T[i], T[1], n0, t_array, rho);
    log_lik[i] = poisson_lpmf(N[i] | yrep[i]); // Log-likelihood calculation
  }
}
