data {
  int<lower=1> S; // Number of steps

  array[S] int<lower=0> N; // observations
  array[S] real T;         // observations
}

parameters {
  real<lower=0> rho_r;        // Parameter rho_r (rate for recoverN)
  real<lower=0> rho_s;        // Parameter rho_s (rate for signal decaN)
  real t0_r;                   // Parameter t_r (time shift)
  real<lower=0, upper=1> f_s;          // Parameter f_s (scaling for the first term)
}

transformed parameters {
  real t_end;
  real ns;
  ns = f_s * N[1];
  t_end = - (1 / rho_s) * log(1 / ns) + T[1];

}

model {
  vector[S] mu;               // Expected values for N given x

  // Define the expected value based on the given equation
  for (i in 1:S) {
    if (T[i] >= t0_r) {
      mu[i] = f_s * N[1] * exp(-rho_s * T[i]) + exp(rho_r * (T[i] - t0_r));
    } else {
      mu[i] = f_s * N[1] * exp(-rho_s * T[i]);
    }

  }

  // Priors
  // n0 ~ normal(N[1], 0.1 * N[1]);
  rho_r ~ normal(0, 1);       // Prior for rho_r
  rho_s ~ normal(0, 1);       // Prior for rho_s
  t0_r ~ normal(0, 1);         // Prior for t_r
  f_s ~ beta(2,2);         // Prior for f_s

  // Likelihood (assuming normallN distributed noise)
  N ~ poisson(mu);
}


