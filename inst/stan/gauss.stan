
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

  int<lower = -1, upper = 1> k; // type of model
  real<lower=0> prior_par;
}

parameters {
  real<lower=0> lambda[k == 1 ? S : 1];
  real<lower=0> mu[k == -1 ? S : 1];
}

model {
  if (k == 0) {
    // weak priors
    lambda[1] ~ inv_gamma(prior_par,prior_par);
    mu[1] ~ inv_gamma(prior_par,prior_par);

    // strong prior for lambda
    // weak prior for mu
    //lambda[1] ~ inv_gamma(2498,12485);
    //mu[1] ~ inv_gamma(23,110);

    for (i in 1:S) {
      if (i == 1) {
        N[i] ~ normal(mean_t(n0, lambda[1], mu[1], T[1]), sigma_t(n0, lambda[1], mu[1], T[1]));
      } else {
        N[i] ~ normal(mean_t(N[i-1], lambda[1], mu[1], T[i] - T[i-1]), sigma_t(N[i-1], lambda[1], mu[1], T[i] - T[i-1]));
      }
    }
  } else if (k == 1) {
    mu[1] ~ inv_gamma(prior_par,prior_par);

    // strong prior
    // mu[1] ~ inv_gamma(2498,12485);

    for (i in 1:S) {
      lambda[i] ~ inv_gamma(prior_par,prior_par);
      // lambda[i] ~ inv_gamma(23, 110); // weak prior

      if (i == 1) {
        N[i] ~ normal(mean_t(n0, lambda[i], mu[1], T[1]), sigma_t(n0, lambda[i], mu[1], T[1]));
      } else {
        N[i] ~ normal(mean_t(N[i-1], lambda[i], mu[1], T[i] - T[i-1]), sigma_t(N[i-1], lambda[i], mu[1], T[i] - T[i-1]));
      }
    }
  } else {
    lambda[1] ~ inv_gamma(prior_par,prior_par);

    // strong prior
    // lambda[1] ~ inv_gamma(2498,12485);

    for (i in 1:S) {
      mu[i] ~ inv_gamma(prior_par,prior_par);

      // weak prior
      // mu[i] ~ inv_gamma(23, 110);
      if (i == 1) {
        N[i] ~ normal(mean_t(n0, lambda[1], mu[i], T[1]), sigma_t(n0, lambda[1], mu[i], T[1]));
      } else {
        N[i] ~ normal(mean_t(N[i-1], lambda[1], mu[i], T[i] - T[i-1]), sigma_t(N[i-1], lambda[1], mu[i], T[i] - T[i-1]));
      }
    }
  }
}

generated quantities {
  real w[k == 0 ? 1 : S];
  real N_rep[S];

  if (k == 0) {
    w[1] = lambda[1] - mu[1];
    for (i in 1:S) {
      if (i == 1) {
        N_rep[i] = normal_rng(mean_t(n0, lambda[1], mu[1], T[1]), sigma_t(n0, lambda[1], mu[1], T[1]));
      } else {
        N_rep[i] = normal_rng(mean_t(N[i-1], lambda[1], mu[1], T[i] - T[i-1]), sigma_t(N[i-1], lambda[1], mu[1], T[i] - T[i-1]));
      }
    }
  } else if (k == -1) {
    for (i in 1:S) {
      w[i] = lambda[1] - mu[i];

      if (i == 1) {
        N_rep[i] = normal_rng(mean_t(n0, lambda[1], mu[i], T[1]), sigma_t(n0, lambda[1], mu[i], T[1]));
      } else {
        N_rep[i] = normal_rng(mean_t(N[i-1], lambda[1], mu[i], T[i] - T[i-1]), sigma_t(N[i-1], lambda[1], mu[i], T[i] - T[i-1]));
      }
    }
  } else {
    for (i in 1:S) {
      w[i] = lambda[i] - mu[1];
      if (i == 1) {
        N_rep[i] = normal_rng(mean_t(n0, lambda[i], mu[1], T[1]), sigma_t(n0, lambda[i], mu[1], T[1]));
      } else {
        N_rep[i] = normal_rng(mean_t(N[i-1], lambda[i], mu[1], T[i] - T[i-1]), sigma_t(N[i-1], lambda[i], mu[1], T[i] - T[i-1]));
      }
    }
  }
}
