
functions {
  real mean_exp (real t, real t0, real[] t_array, real[] rho_array) {
    real res;
    int n_t = num_elements(t_array);
    int n_rho = num_elements(rho_array);

    if (n_t == 0) return exp(rho_array[1] * (t - t0));

    if (t <= t_array[1]) return exp(rho_array[1] * (t - t0));

    res = exp(rho_array[1] * (t_array[1] - t0));

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

   real logistic_growth(real t, real n0, real rho, real K) {
    real num = K * n0;
    real den = n0 + (K - n0) * exp(-rho * t);
    return(num/den);
  }

  real mean_log (real t, real t0, real[] t_array, real[] rho_array, real K) {
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

  real<lower=0> prior_K;
}

parameters {
  real omega;
  real rho_exp[G];
  real rho_log[G];
  real<upper=T[1]> t0_exp;
  real<upper=T[1]> t0_log;
  real<lower=prior_K> K; // carrying capacity
}

model {
  target += beta_lpdf(omega | 5, 5);
  target += normal_lpdf(K | prior_K, prior_K); // sample the carrying capacity

  for (i in 1:G) {
    target += normal_lpdf(rho_exp[i] | 0, 1);
    target += normal_lpdf(rho_log[i] | 0, 1);
  }

  target += normal_lpdf(t0_exp | T[1], 100);
  target += normal_lpdf(t0_log | T[1], 100);

  for (i in 1:S) {
    target += log_mix(
      omega,
      poisson_lpmf(N[i] | mean_exp(T[i], t0_exp, t_array, rho_exp)),
      poisson_lpmf(N[i] | mean_log(T[i], t0_log, t_array, rho_log, K))
    );
  }
}
