
functions {
  real mean_t (real t, real t0, real t0_r, real ns, real rho_s, real rho_r) {
    real s_pop;
    real r_pop;

    s_pop = ns * exp(-rho_s * (t - t0));
    if (t >= t0_r) {
      r_pop = 1 * exp(rho_r * (t - t0_r));
    } else {
      r_pop = 0;
    }

    return s_pop + r_pop;
  }
}

data {
  int<lower=1> S; // Number of steps

  array[S] int<lower=0> N; // observations
  array[S] real T;         // observations
}

parameters {
  real<lower=0> rho_s;
  real<lower=0> rho_r;

  real<lower=1, upper=N[1] * 1.1> ns;
  real t0_r;
}

transformed parameters {
  real t_end;
  t_end = - (1 / rho_s) * log(1 / ns) + T[1];
}

model {
  target += inv_gamma_lpdf(rho_r | 1, 1);
  target += inv_gamma_lpdf(rho_s | 1, 1);

  target += normal_lpdf(ns | (1 + N[1]) / 2.0, N[1]);
  target += normal_lpdf(t0_r | T[1], 100);

  for (i in 1:S) {
    target += poisson_lpmf(N[i] | mean_t(T[i], T[1], t0_r, ns, rho_s, rho_r));
  }
}
