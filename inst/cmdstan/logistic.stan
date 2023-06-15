
functions {
  real logistic_growth(real t, real n0, real rho, real K) {
    real num = K * n0;
    real den = n0 + (K - n0) * exp(-rho * t);
    return(num/den);
  }

  real mean_t (real t, real t0, real[] t_array, real[] rho_array, real K) {
    real current_n0 = 1;
    real dt;
    int n_t = num_elements(t_array);
    int n_rho = num_elements(rho_array);

    if (n_t == 0) {
      dt = t- t0;
      return(logistic_growth(dt, current_n0, rho_array[1], K));
    }

    if (t <= t_array[1]) {
      dt = t- t0;
      return(logistic_growth(dt, current_n0, rho_array[1], K));
    }

    dt = t_array[1] - t0;
    current_n0 = logistic_growth(dt, current_n0, rho_array[1], K);
    for (i in 2:n_t) {
      if (t <= t_array[i]) {
        dt = t - t_array[i-1];
        return(logistic_growth(dt, current_n0, rho_array[i], K));
      } else {
        dt = t_array[i] - t_array[i-1];
        current_n0 = logistic_growth(dt, current_n0, rho_array[i], K);
      }
    }

    dt = t - t_array[n_t];
    return(logistic_growth(dt, current_n0, rho_array[n_rho], K));
  }
}

data {
  int<lower=1> S; // Number of steps
  int<lower=1> G; // Number of wondows

  int <lower=0> N[S];      // observations
  real T[S];      // observations

  real t_array[G - 1];

  real<upper=T[1]> t0_lower_bound;
  real<lower=0> prior_K;
}

parameters {
  real rho[G];
  real<upper=T[1]> t0; // t0 will be smaller or equal than the first time point we have
  // real<lower=prior_K, upper=100 * prior_K> K; // carrying capacity
  real<lower=prior_K> K; // carrying capacity
}

model {
  // target += uniform_lpdf(K | prior_K, 100 * prior_K); // sample the carrying capacity
  target += normal_lpdf(K | prior_K, prior_K); // sample the carrying capacity
  for (i in 1:G) {
    target += normal_lpdf(rho[i] | 0, 1); // sample the growth rates
  }

  // target += normal_lpdf(t0 | - (1 / rho[1]) * log((K - 1 - N[1]) / (N[1] * (K - 1))), 1);// sample t0
  target += uniform_lpdf(t0 | t0_lower_bound, T[1]);

  for (i in 1:S) {
    target += poisson_lpmf(N[i] | mean_t(T[i], t0, t_array, rho, K));
  }
}

// generated quantities {
//   int <lower=0> N_rep[S];
//   real log_lik[S];
//
//   for (i in 1:S) {
//     log_lik[i] = poisson_lpmf(N[i] | mean_t(T[i], t0, t_array, rho, K));
//     N_rep[i] = poisson_rng(mean_t(T[i], t0, t_array, rho, K));
//   }
// }

