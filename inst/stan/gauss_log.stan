
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

  int<lower = -1, upper = 1> k; // type of model
  real<lower=0> prior_par;
  real<lower=0> max_K;
}

parameters {
  real<lower=0> lambda[1];
  real<lower=0> mu[k == 1 ? S : 1];
  real<lower=max_K> K; // carrying capacity
}

model {
  K ~ normal(max_K, max_K * .5);
  lambda[1] ~ inv_gamma(prior_par,prior_par);

  if (k == 0) {
    // weak prior
    mu[1] ~ inv_gamma(prior_par,prior_par);

    for (i in 1:S) {
      N[i] ~ normal(mean_t(n0, lambda[1], mu[1], K, T[i]), sigma_t(n0, lambda[1], mu[1], K, T[i]));
    }

  } else {

    for (i in 1:S) {
      mu[i] ~ inv_gamma(prior_par, prior_par);

      if (i == 1) {
        N[i] ~ normal(mean_t(n0, lambda[1], mu[i], K, T[1]), sigma_t(n0, lambda[1], mu[i], K, T[1]));
      } else {
        N[i] ~ normal(mean_t(N[i-1], lambda[1], mu[i], K, T[i] - T[i-1]), sigma_t(N[i-1], lambda[1], mu[i], K, T[i] - T[i-1]));
      }
    }
  }
}

generated quantities {
  real ro[k == 0 ? 1 : S];
  real N_rep[S];

  if (k == 0) {
    ro[1] = lambda[1] - mu[1];
    for (i in 1:S) {
      N_rep[i] = normal_rng(mean_t(n0, lambda[1], mu[1], K, T[i]), sigma_t(n0, lambda[1], mu[1], K, T[i]));
    }
  } else {
    for (i in 1:S) {
      ro[i] = lambda[1] - mu[i];
      if (i == 1) {
        N_rep[i] = normal_rng(mean_t(n0, lambda[1], mu[i], K, T[1]), sigma_t(n0, lambda[1], mu[i], K, T[1]));
      } else {
        N_rep[i] = normal_rng(mean_t(N[i-1], lambda[1], mu[i], K, T[i] - T[i-1]), sigma_t(N[i-1], lambda[1], mu[i], K, T[i] - T[i-1]));
      }
    }
  }
}


