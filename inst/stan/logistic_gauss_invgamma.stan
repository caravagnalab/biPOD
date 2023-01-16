
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
  int<lower=1> S; // Number of steps
  real<lower=0> n0;
  real<lower=0> t0;

  real N[S];   // observations
  real T[S];  // observations

  real<lower=0> a;
  real<lower=0> b;
  real<lower=0> prior_K;
}

parameters {
  real<lower=0> lambda;
  real<lower=0> mu;
  real<lower=0> K; // carrying capacity
}

transformed parameters {
  real ro;
  ro = lambda - mu;
}

model {
  K ~ normal(prior_K, prior_K * .1);
  lambda ~ inv_gamma(a, b);
  mu ~ inv_gamma(a, b);

  for (i in 1:S) {
    N[i] ~ normal(mean_t(n0, lambda, mu, K, T[i] - t0), sigma_t(n0, lambda, mu, K, T[i] - t0));
  }
}

generated quantities {
  real N_rep[S];

  for (i in 1:S) {
    N_rep[i] = normal_rng(mean_t(n0, lambda, mu, K, T[i] - t0), sigma_t(n0, lambda, mu, K, T[i] - t0));
  }
}


