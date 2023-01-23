
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
  int<lower=1> S; // Number of steps
  real<lower=0> n0;
  real<lower=0> t0;

  real N[S];   // observations
  real T[S];  // observations

  real<lower=0> b;
  real<lower=0> a;
}

parameters {
  real <lower=0> lambda;
  real <lower=0> mu;
}

transformed parameters {
  real ro;
  ro = lambda - mu;
}

model {
  target += inv_gamma_lpdf(lambda | a, b);
  target += inv_gamma_lpdf(mu | a, b);

  for (i in 1:S) {
    target += normal_lpdf(N[i] | mean_t(n0, lambda, mu, T[i] - t0), sigma_t(n0, lambda, mu, T[i] - t0));
  }
}

generated quantities {
  real N_rep[S];
  real log_lik[S];

  for (i in 1:S) {
    log_lik[i] = normal_lpdf(N[i] | mean_t(n0, lambda, mu, T[i] - t0), sigma_t(n0, lambda, mu, T[i] - t0));
    N_rep[i] = normal_rng(mean_t(n0, lambda, mu, T[i] - t0), sigma_t(n0, lambda, mu, T[i] - t0));
  }
}
