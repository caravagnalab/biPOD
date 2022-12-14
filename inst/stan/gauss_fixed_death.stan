
functions {
  real mean_t (real n0, real lambda, real mu, real t) {
    real m;
    m = n0 * exp(t * (lambda - mu));
    return m;
  }

  real sigma_t (real n0, real lambda, real mu, real t) {
    real s;
    real w = lambda - mu;

    if (lambda == mu) {
      s = n0 * 2 * lambda * t;
    } else {
      s = n0 * ((lambda + mu) / w) * exp(w*t) * (exp(w*t) - 1);
    }
    return sqrt(s);
  }
}

data {
  int<lower=2> S; // Number of steps
  real<lower=0> n0;

  real N[S];   // observations
  real T[S];  // observations
}

parameters {
  real<lower=0> lambda[S]; // birth rates are variable
  real<lower=0> mu;        // death rate is fixed
}

model {
  mu ~ inv_gamma(2,2);

  for (i in 1:S) {
    lambda[i] ~ inv_gamma(2,2);
    if (i == 1) {
      N[i] ~ normal(mean_t(n0, lambda[i], mu, T[1]), sigma_t(n0, lambda[i], mu, T[1]));
    } else {
      N[i] ~ normal(mean_t(N[i-1], lambda[i], mu, T[i] - T[i-1]), sigma_t(N[i-1], lambda[i], mu, T[i] - T[i-1]));
    }
  }
}

generated quantities {
  real N_rep[S];
  for (i in 1:S) {
    if (i == 1) {
      N_rep[i] = normal_rng(mean_t(n0, lambda[i], mu, T[1]), sigma_t(n0, lambda[i], mu, T[1]));
    } else {
      N_rep[i] = normal_rng(mean_t(N[i-1], lambda[i], mu, T[i] - T[i-1]), sigma_t(N[i-1], lambda[i], mu, T[i] - T[i-1]));
    }
  }
}
