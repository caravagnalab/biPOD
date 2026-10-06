functions {
  // Populations on the log scale: exp() of the exponents overflows on long
  // follow-up. An absent population keeps the floor of 1e-9 (a mean of
  // exactly 0 would make any positive count impossible and leave hard walls
  // in the posterior).
  real log_ns(real t, real rho_s, real t_end) {
    return t < t_end ? -rho_s * (t - t_end) : log(1e-9);
  }

  real log_nr(real t, real rho_r, real t0_r) {
    return t >= t0_r ? rho_r * (t - t0_r) : log(1e-9);
  }
}

data {
  int<lower=1> S; // Number of steps
  array[S] int<lower=0> N; // observations
  array[S] real T;         // observations
  int<lower=0,upper=1> prior_only;
}

parameters {
  real<lower=0> rho_r;        // Parameter rho_r (rate for recover)
  real<upper=T[1]> t0_r;                   // Parameter t_r (time shift)
}

model {
  vector[S] log_mu;           // log expected values for N

  // Priors
  rho_r ~ normal(0, 1);       // Prior for rho_r
  t0_r ~ normal(T[1], T[S] - T[1]);         // Prior for t_r

  for (i in 1:S) {
    log_mu[i] = log_nr(T[i], rho_r, t0_r);
  }
  N ~ poisson_log(log_mu);
}

generated quantities {
  vector[S] log_lik;             // Log-likelihood for each observation
  vector[S] yrep;               // Expected values for N given x
  vector[S] ns;
  vector[S] nr;

  for (i in 1:S) {
    nr[i] = exp(log_nr(T[i], rho_r, t0_r));
    ns[i] = 0.0;
    yrep[i] = nr[i];
    log_lik[i] = poisson_log_lpmf(N[i] | log_nr(T[i], rho_r, t0_r));
  }
}
