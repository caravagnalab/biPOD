
functions {
  real birthDeathLike_log (int m, int n, real delta_t, real lambda, real mu) {
    real prob;
    real lprob;
    real a_t;
    real b_t;
    int upper_v;
    int lower_v;

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
  real<lower=0> lambda; // birth rate is fixed
  real<lower=0> mu;     // death rate is fixed
}

model {
  lambda ~ inv_gamma(2,2);
  mu ~ inv_gamma(2,2);

  for (i in 1:S) {
    if (i == 1) {
      N[i] ~ birthDeathLike(n0, T[i], lambda, mu);
    } else {
      N[i] ~ birthDeathLike(N[i-1], T[i] - T[i-1], lambda, mu);
    }
  }
}

