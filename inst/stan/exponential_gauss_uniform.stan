
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
  real T[S];   // time-points

  real<lower=0> b;
  real<lower=0> a;
  real<lower=0> g;
}

parameters {
  real <lower=a, upper = b> lambda;
  real <lower=a, upper = b> mu;
}

transformed parameters {
  real ro;
  ro = lambda - mu;
}

model {
  // Use uniform priors sampled on previous centered value
  target += uniform_lpdf(lambda | lambda - (b-a)/g, lambda + (b-a)/g);
  target += uniform_lpdf(mu | mu - (b-a)/g, mu + (b-a)/g);

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
