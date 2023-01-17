
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

    psi = (n0 * exp((lambda-mu) * dt) + K - n0) / (K * exp((lambda - mu) * dt));
    phi = (K - n0) * (mu + lambda - mu) * (1 - exp(-dt * (lambda - mu))) / ((lambda - mu) * K) + mu * n0 * dt / K;
    s = n0 * (psi + 2 * phi - 1) / psi^2;

    return sqrt(s);
  }
}

data {
  int<lower=1> S; // Number of steps
  real<lower=0> n0;
  real<lower=0> t0;

  int N[S];   // observations
  real T[S];  // observations

  real<lower=0> a;
  real<lower=0> b;
  real<lower=0> g;
  real<lower=0> prior_K;
}

parameters {
  real<lower=a, upper=b> lambda;
  real<lower=a, upper=b> mu;
  real<lower=prior_K * .95> K; // carrying capacity
}

transformed parameters {
  real ro;
  ro = lambda - mu;
}

model {
  K ~ normal(prior_K, prior_K * .01);
  lambda ~ uniform(lambda - (b-a)/g, lambda + (b-a)/g);
  mu ~ uniform(mu - (b-a)/g, mu + (b-a)/g);

  for (i in 1:S) {
    // N[i] ~ normal(mean_t(n0, lambda, mu, K, T[i] - t0), sigma_t(n0, lambda, mu, K, T[i] - t0));
    N[i] ~ poisson(mean_t(n0, lambda, mu, K, T[i] - t0));
  }
}

generated quantities {
  int N_rep[S];

  for (i in 1:S) {
    // N_rep[i] = normal_rng(mean_t(n0, lambda, mu, K, T[i] - t0), sigma_t(n0, lambda, mu, K, T[i] - t0));
    N_rep[i] = poisson_rng(mean_t(n0, lambda, mu, K, T[i] - t0));
  }
}
