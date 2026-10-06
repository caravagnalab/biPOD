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
  real<lower=0> rho_r;        // Parameter rho_r (rate for recoverN)
  real<lower=0> rho_s;        // Parameter rho_s (rate for signal decaN)
  real<lower=T[1]> t0_r;                   // Parameter t_r (time shift)
  real<lower=T[1]> t_end;
}

model {
  vector[S] log_mu;           // log expected values for N

  // Priors
  rho_r ~ normal(0, 1);       // Prior for rho_r
  rho_s ~ normal(0, 1);       // Prior for rho_s
  t0_r ~ normal(T[1], T[S] - T[1]);         // Prior for t_r
  t_end ~ normal(T[S], T[S] - T[1]);

  for (i in 1:S) {
    log_mu[i] = log_sum_exp(log_ns(T[i], rho_s, t_end), log_nr(T[i], rho_r, t0_r));
  }
  N ~ poisson_log(log_mu);
}

generated quantities {
  vector[S] log_lik;             // Log-likelihood for each observation
  vector[S] ns;
  vector[S] nr;
  vector[S] yrep;               // Expected values for N given x

  for (i in 1:S) {
    real log_mu = log_sum_exp(log_ns(T[i], rho_s, t_end), log_nr(T[i], rho_r, t0_r));
    ns[i] = exp(log_ns(T[i], rho_s, t_end));
    nr[i] = exp(log_nr(T[i], rho_r, t0_r));
    yrep[i] = exp(log_mu);
    log_lik[i] = poisson_log_lpmf(N[i] | log_mu);
  }
}
