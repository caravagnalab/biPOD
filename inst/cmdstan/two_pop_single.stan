data {
  int<lower=1> S; // Number of steps
  array[S] int<lower=0> N; // observations
  array[S] real T;         // observations
}

parameters {
  real<lower=0> rho_r;        // Parameter rho_r (rate for recover)
  real<upper=T[1]> t0_r;                   // Parameter t_r (time shift)
}

model {
  vector[S] mu;               // Expected values for N given x
  vector[S] ns;
  vector[S] nr;

  // Define the expected value based on the given equation
  for (i in 1:S) {
    nr[i] = exp(rho_r * (T[i] - t0_r));
    ns[i] = 0.0;
    mu[i] =  nr[i] + ns[i];
  }

  // Priors
  rho_r ~ normal(0, 1);       // Prior for rho_r
  t0_r ~ normal(T[1], 5);         // Prior for t_r

  // Likelihood (assuming normallN distributed noise)
  N ~ poisson(mu);
}

generated quantities {
  vector[S] log_lik;             // Log-likelihood for each observation
  vector[S] yrep;               // Expected values for N given x
  vector[S] ns;
  vector[S] nr;

  for (i in 1:S) {
    nr[i] = exp(rho_r * (T[i] - t0_r));
    ns[i] = 0.0;
    yrep[i] =  nr[i] + ns[i];
    log_lik[i] = poisson_lpmf(N[i] | yrep[i]); // Log-likelihood calculation
  }
}
