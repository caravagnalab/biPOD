data {
  int<lower=1> S; // Number of steps
  array[S] int<lower=0> N; // observations
  array[S] real T;         // observations
  int<lower=0,upper=1> prior_only;
}

parameters {
  real<lower=0> rho_s;        // Parameter rho_s (rate for decay)
  real<lower=T[1]> t_end;     // Parameter t_end (time shift)
  real<lower=0> phi;           // NB overdispersion
}

model {
  vector[S] mu;               // Expected values for N given x
  vector[S] ns;
  vector[S] nr;

  for (i in 1:S) {
    if (T[i] >= t_end) {
      ns[i] = 1e-9;
    } else {
      ns[i] = exp(-rho_s * (T[i] - t_end));
    }
    nr[i] = 0.0;
    mu[i] = ns[i] + nr[i];
  }

  // Priors
  rho_s ~ normal(0, 1);       // Prior for rho_s
  t_end ~ normal(T[S], T[S] - T[1]);    // Prior for t_end
  phi ~ gamma(2, 0.1);
  N ~ neg_binomial_2(mu, phi);
}

generated quantities {
  vector[S] log_lik;             // Log-likelihood for each observation
  vector[S] yrep;                // Expected values for N given x
  vector[S] ns;
  vector[S] nr;

  for (i in 1:S) {
    if (T[i] >= t_end) {
      ns[i] = 1e-9;
    } else {
      ns[i] = exp(-rho_s * (T[i] - t_end));
    }
    nr[i] = 0.0;
    yrep[i] = ns[i] + nr[i];
    log_lik[i] = neg_binomial_2_lpmf(N[i] | yrep[i], phi); // Log-likelihood calculation
  }
}
