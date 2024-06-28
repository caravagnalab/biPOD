functions {
  real ind(real x, real y) {
    if (x >= y) {
      return 1.0;
    } else {
      return 0.0;
    }
  }

  real expected_mean(real x, real q, vector s, vector b) {
    int G = num_elements(s);
    real res = q + x * s[1];
    for (g in 2:G) {
      res = res + (x - b[g-1]) * s[g] * ind(x, b[g-1]);
    }
    return exp(res);
  }
}

data {
  int<lower=1> S; // Number of steps
  int<lower=0> G; // Number of windows

  array[S] real N; // observations
  array[S] real T;         // observations

  vector[G] b;
}

parameters {
  real q; // intercept
  vector[G+1] s; // slopes
  real<lower=0> sigma;
}

model {
  target += normal_lpdf(q | 0, 1);
  target += inv_gamma_lpdf(sigma | .001, .001);

  for (g in 1:(G+1)) {
    target += normal_lpdf(s[g] | 0, 1);
  }

  for (i in 1:S) {
    target += normal_lpdf(N[i] | expected_mean(T[i], q, s, b), sigma);
  }
}

generated quantities {
  vector[S] log_lik;
  for (i in 1:S) {
    log_lik[i] = normal_lpdf(N[i] | expected_mean(T[i], q, s, b), sigma);
  }
}
