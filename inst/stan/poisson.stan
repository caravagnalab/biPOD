
data {
  int<lower=2> S; // Number of steps
  int<lower=0> n0;

  int N[S];   // observations
  real T[S];  // observations

  int<lower = -1, upper = 1> k; // type of model
  real<lower=0> prior_par;
}

parameters {
  real ro[k == 0 ? 1 : S];
}

model {
  if (k == 0) {
    ro[1] ~ normal(0, prior_par);

    for (i in 1:S) {
      if (i == 1) {
        N[i] ~ poisson(n0 * exp(ro[1] * T[1]));
      } else {
        N[i] ~ poisson(N[i-1] * exp(ro[1] * (T[i] - T[i-1])));
      }
    }
  } else {

    for (i in 1:S) {
      // priors
      ro[i] ~ normal(0, prior_par);

      if (i == 1) {
        N[i] ~ poisson(n0 * exp(ro[i] * T[1]));
      } else {
        N[i] ~ poisson(N[i-1] * exp(ro[i] * (T[i] - T[i-1])));
      }
    }
  }
}

generated quantities {
  int<lower=0> N_rep[S];

  if (k == 0) {
    for (i in 1:S) {
      if (i == 1) {
        N_rep[i] = poisson_rng(n0 * exp(ro[1] * T[1]));
      } else {
        N_rep[i] = poisson_rng(N[i-1] * exp(ro[1] * (T[i] - T[i-1])));
      }
    }
  } else {
    for (i in 1:S) {
      if (i == 1) {
        N_rep[i] = poisson_rng(n0 * exp(ro[i] * T[1]));
      } else {
        N_rep[i] = poisson_rng(N[i-1] * exp(ro[i] * (T[i] - T[i-1])));
      }
    }
  }
}
