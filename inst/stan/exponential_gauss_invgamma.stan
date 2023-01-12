
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

  int<lower=0, upper=1> k; // type of model
  real<lower=0> b;
  real<lower=0> a;
}

parameters {
  real <lower=0> lambda;
  real <lower=0> mu[k == 0 ? 1 : S];
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
  lambda ~ inv_gamma(a,b);
  if (k == 0) {
    mu[1] ~ inv_gamma(a,b);
  } else {
    for (i in 1:S) {
      mu[i] ~ inv_gamma(a,b);
    }
  }

  if (k == 0) {
    for (i in 1:S) {
      N[i] ~ normal(mean_t(n0, lambda, mu[1], T[i]), sigma_t(n0, lambda, mu[1], T[i]));
    }
  } else {

    for (i in 1:S) {
      if (i == 1) {
        N[i] ~ normal(mean_t(n0, lambda, mu[1], T[1]), sigma_t(n0, lambda, mu[1], T[1]));
      } else {
        N[i] ~ normal(mean_t(N[i-1], lambda, mu[i], T[i] - T[i-1]), sigma_t(N[i-1], lambda, mu[i], T[i] - T[i-1]));
      }
    }
  }
}

generated quantities {
  real N_rep[S];

  if (k == 0) {
    for (i in 1:S) {
      if (i == 1) {
        N_rep[i] = normal_rng(mean_t(n0, lambda, mu[1], T[1]), sigma_t(n0, lambda, mu[1], T[1]));
      } else {
        N_rep[i] = normal_rng(mean_t(N[i-1], lambda, mu[1], T[i] - T[i-1]), sigma_t(N[i-1], lambda, mu[1], T[i] - T[i-1]));
      }
    }
  } else {
    for (i in 1:S) {
      if (i == 1) {
        N_rep[i] = normal_rng(mean_t(n0, lambda, mu[i], T[1]), sigma_t(n0, lambda, mu[i], T[1]));
      } else {
        N_rep[i] = normal_rng(mean_t(N[i-1], lambda, mu[i], T[i] - T[i-1]), sigma_t(N[i-1], lambda, mu[i], T[i] - T[i-1]));
      }
    }
  }
}
