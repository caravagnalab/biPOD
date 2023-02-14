
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
    real w = lambda - mu;

    phi = (K - n0) * (lambda) / (K * (lambda - mu)) * (1 - exp(-(lambda - mu) * dt)) + mu * n0 * dt / K;
    psi = (n0 * exp((lambda - mu) * dt) + K - n0) / (K * exp((lambda - mu) * dt));
    s = n0 * (psi + 2 * phi - 1) / (psi^2);

    return sqrt(s);

    // if (lambda == mu) {
    //   s = n0 * 2 * lambda * dt;
    // } else {
    //   s = n0 * ((lambda + mu) / w) * exp(w * dt) * (exp(w * dt) - 1);
    // }
    // return sqrt(s);
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
  real<lower=0> g;
  real<lower=0> prior_K;
}

parameters {
  real <lower=a, upper = b> lambda;
  real <lower=a, upper = b> mu;
  real<lower=prior_K * .95> K; // carrying capacity
}

transformed parameters {
  real ro;
  ro = lambda - mu;
}

model {
  target += normal_lpdf(K | prior_K, prior_K * .1);
  target += uniform_lpdf(lambda | lambda - (b-a)/g, lambda + (b-a)/g);
  target += uniform_lpdf(mu | mu - (b-a)/g, mu + (b-a)/g);

  for (i in 1:S) {
    target += normal_lpdf(N[i] | mean_t(n0, lambda, mu, K, T[i] - t0), sigma_t(n0, lambda, mu, K, T[i] - t0));
  }
}

generated quantities {
  real N_rep[S];
  real log_lik[S];

  for (i in 1:S) {
    log_lik[i] = normal_lpdf(N[i] | mean_t(n0, lambda, mu, K, T[i] - t0), sigma_t(n0, lambda, mu, K, T[i] - t0));
    N_rep[i] = normal_rng(mean_t(n0, lambda, mu, K, T[i] - t0), sigma_t(n0, lambda, mu, K, T[i] - t0));
  }
}
