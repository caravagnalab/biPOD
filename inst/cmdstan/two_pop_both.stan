data {
  int<lower=1> S; // Number of steps
  array[S] int<lower=0> N; // observations
  array[S] real T;         // observations
}

parameters {
  real<lower=0> rho_r;        // Parameter rho_r (rate for recoverN)
  real<lower=0> rho_s;        // Parameter rho_s (rate for signal decaN)
  real t0_r;                   // Parameter t_r (time shift)
  real<lower=T[1]> t_end;
  // real<lower=0, upper=1.5> f_s;
}


model {
  vector[S] mu;               // Expected values for N given x
  vector[S] ns;
  vector[S] nr;

  // Priors
  // n0 ~ normal(N[1], 0.1 * N[1]);
  rho_r ~ normal(0, 1);       // Prior for rho_r
  rho_s ~ normal(0, 1);       // Prior for rho_s
  t0_r ~ normal(T[1], 5);         // Prior for t_r
  t_end ~ normal(T[S], 5);

  for (i in 1:S) {
    if (T[i] >= t_end) {
      ns[i] = 1e-9;
    } else {
      ns[i] = exp(-rho_s * (T[i] - t_end));
    }

    if (T[i] >= t0_r) {
      nr[i] = exp(rho_r * (T[i] - t0_r));
    } else {
      nr[i] = 1e-9;
    }

    mu[i] = nr[i] + ns[i];
  }

  // Likelihood (assuming normallN distributed noise)
  N ~ poisson(mu);
}

generated quantities {
  vector[S] log_lik;             // Log-likelihood for each observation
  vector[S] ns;
  vector[S] nr;
  vector[S] yrep;               // Expected values for N given x

  for (i in 1:S) {
    if (T[i] >= t_end) {
      ns[i] = 1e-9;
    } else {
      ns[i] = exp(-rho_s * (T[i] - t_end));
    }

    if (T[i] >= t0_r) {
      nr[i] = exp(rho_r * (T[i] - t0_r));
    } else {
      nr[i] = 1e-9;
    }
    yrep[i] = nr[i] + ns[i];
    log_lik[i] = poisson_lpmf(N[i] | yrep[i]); // Log-likelihood calculation
  }
}
