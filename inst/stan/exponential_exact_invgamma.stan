
functions {
  real birthDeathLike_log (int m, int n, real delta_t, real lambda, real mu) {
    real prob;
    real lprob;
    real alpha_t;
    real beta_t;
    int upper_sum_value;

    if (n <= m) {
      upper_sum_value = n;
    } else {
      upper_sum_value = m;
    }

    alpha_t = (mu * (exp(delta_t*(lambda - mu)) - 1)) / (lambda * exp((lambda - mu) * delta_t) - mu);
    beta_t = lambda * alpha_t / mu;

    prob = 0;

    for (j in 1:upper_sum_value) {
      prob = prob + exp(lchoose(n,j)) * (1 - alpha_t)^j * alpha_t^(n-j) * exp(lchoose(m-1, m-j)) * (1 - beta_t)^j * beta_t^(m-j);
    }

    lprob = log(prob);
    return lprob;
  }

  int birthDeathLike_rng (int n, real delta_t, real lambda, real mu) {
    real prob;
    real lprob;
    real a_t;
    real b_t;
    int upper_v;
    int lower_v;
    real x;
    real cdf;
    int m;

    x = uniform_rng(0, 1);
    cdf = 0;
    m = -1;

    while (cdf < x) {
      m += 1;

      upper_v = n - 1;
      if (n - m >= 0) {
        lower_v = n - m;
      } else {
        lower_v = 0;
      }

      a_t = (mu * (exp(delta_t*(lambda - mu)) - 1)) / (lambda * exp((lambda - mu) * delta_t) - mu);
      b_t = lambda * a_t / mu;

      prob = 0;

      for (j in lower_v:upper_v) {
        prob = prob + exp(lchoose(n, j)) * exp(lchoose(m - 1, n - j -1)) * a_t^j * ((1-a_t)*(1-b_t))^(n-j) * b_t^(m-n+j);
      }

      cdf += prob;
    }

    return m;
  }
}

data {
  int<lower=2> S; // Number of steps
  int<lower=0> n0;

  int N[S];   // observations
  real T[S];  // observations

  int<lower=0, upper=1> k; // type of model
  real<lower=0> a;
  real<lower=0> b;
}

parameters {
  real<lower=0> lambda;
  real<lower=0> mu[k == 1 ? S : 1];
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
      N[i] ~ birthDeathLike(n0, T[i], lambda, mu[1]);
    }
  } else {
    for (i in 1:S) {
      if (i == 1) {
        N[i] ~ birthDeathLike(n0, T[i], lambda, mu[i]);
      } else {
        N[i] ~ birthDeathLike(N[i-1], T[i] - T[i-1], lambda, mu[i]);
      }
    }
  }
}

generated quantities {}
