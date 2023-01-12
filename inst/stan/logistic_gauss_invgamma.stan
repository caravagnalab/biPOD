
functions {
  real mean_t (real n0, real lambda, real mu, real K, real dt) {
    real m;

    m = (K * n0) / (n0 + (K - n0) * exp(- (lambda - mu) * dt));

    return m;
  }

  real sigma_t (real n0, real lambda, real mu, real K, real dt) {
    real s;
    real phi;
    real psi;

    phi = (n0 * exp((lambda - mu) * dt) + K - n0) / (K * exp((lambda - mu) * dt));
    psi = (K - n0) * (lambda) / (K * (lambda - mu)) * (1 - exp(-(lambda - mu) * dt)) + (mu + n0 * dt) / K;
    s = n0 * (phi + 2*psi - 1) / (phi^2);

    return sqrt(s);
  }
}

data {
  int<lower=2> S; // Number of steps
  real<lower=0> n0;

  real N[S];   // observations
  real T[S];  // observations

  int<lower = 0, upper = 1> k; // type of model
  real<lower=0> a;
  real<lower=0> b;
  real<lower=0> prior_K;
}

parameters {
  real<lower=0> lambda;
  real<lower=0> mu[k == 1 ? S : 1];
  real K; // carrying capacity
}

transformed parameters {
  real ro[k == 0 ? 1 : S];
  if (k == 0) {
    ro[1] = lambda - mu[1];
  } else {
    for (i in 1:S) {
      ro[i] = lambda - mu[i];
    }
  }
}

model {
  K ~ normal(prior_K, prior_K * .5);
  lambda ~ inv_gamma(a, b);

  if (k == 0) {
    mu[1] ~ inv_gamma(a, b);
  } else {
    for (i in 1:S) {
      mu[i] ~ inv_gamma(a, b);
    }
  }

  if (k == 0) {
    for (i in 1:S) {
      N[i] ~ normal(mean_t(n0, lambda, mu[1], K, T[i]), sigma_t(n0, lambda, mu[1], K, T[i]));
    }
  } else {
    for (i in 1:S) {
      if (i == 1) {
        N[i] ~ normal(mean_t(n0, lambda, mu[i], K, T[1]), sigma_t(n0, lambda, mu[i], K, T[1]));
      } else {
        N[i] ~ normal(mean_t(N[i-1], lambda, mu[i], K, T[i] - T[i-1]), sigma_t(N[i-1], lambda, mu[i], K, T[i] - T[i-1]));
      }
    }
  }
}

generated quantities {
  real N_rep[S];

  if (k == 0) {
    for (i in 1:S) {
      N_rep[i] = normal_rng(mean_t(n0, lambda, mu[1], K, T[i]), sigma_t(n0, lambda, mu[1], K, T[i]));
    }
  } else {
    for (i in 1:S) {
      if (i == 1) {
        N_rep[i] = normal_rng(mean_t(n0, lambda, mu[i], K, T[1]), sigma_t(n0, lambda, mu[i], K, T[1]));
      } else {
        N_rep[i] = normal_rng(mean_t(N[i-1], lambda, mu[i], K, T[i] - T[i-1]), sigma_t(N[i-1], lambda, mu[i], K, T[i] - T[i-1]));
      }
    }
  }
}


