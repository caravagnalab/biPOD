
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
      if (exp(lchoose(n,j)) != positive_infinity() && exp(lchoose(m-1,m-j)) != positive_infinity()) {
        prob = prob + exp(lchoose(n,j)) * (1 - alpha_t)^j * alpha_t^(n-j) * exp(lchoose(m-1, m-j)) * (1 - beta_t)^j * beta_t^(m-j);
      }
    }

    lprob = log(prob);
    return lprob;
  }
}

data {
  int<lower=2> S; // Number of steps
  int<lower=0> n0;

  int N[S];   // observations
  real T[S];  // observations
}

parameters {
  real<lower=0> lambda[S]; // birth rates are variable
  real<lower=0> mu[S];     // death rates are variable
}

model {
  for (i in 1:S) {
    lambda[i] ~ inv_gamma(2,2);
    mu[i] ~ inv_gamma(2,2);
    if (i == 1) {
      N[i] ~ birthDeathLike(n0, T[i], lambda[i], mu[i]);
    } else {
      N[i] ~ birthDeathLike(N[i-1], T[i] - T[i-1], lambda[i], mu[i]);
    }
  }
}

