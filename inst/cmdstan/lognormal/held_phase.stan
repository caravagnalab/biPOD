// Held-phase model for a continuous marker with a non-tumour background:
// decline, optional held phase, regrowth.
// Noise is Student-t on the log scale (kept under lognormal/ with the other
// continuous-marker models).
//
// Tumour-derived marker T(t) on the log scale:
//   log T(t) = log_T0 + k_d * (min(t, t_n) - t0) + r * max(t - t_e, 0)
// Decline at rate k_d < 0 until t_n; flat between t_n and t_e; regrowth at
// r > 0 after t_e. The decline needs an intervention, so t_n <= t_n_max (its
// end + lag). held = 0 fixes t_e = t_n (continuous regrowth from the nadir,
// model M0); held = 1 frees t_e in [t_n, T[S]] (model M1).
//
// Observed marker = background b + T(t), Student-t noise on the log scale:
//   log y ~ student_t(nu, log(b + T(t)), sigma)

data {
  int<lower=2> S;
  vector[S] t;                    // sampling times (days), sorted
  vector[S] y;                    // marker values
  real<lower=t[1]> t_n_max;
  int<lower=0, upper=1> held;
  vector<lower=0, upper=1>[S] w;  // likelihood weights: 0 leaves a sample out (exact LOO refits)
  real log_b_prior;
  real<lower=0> log_b_sd;         // prior sd of log background (0.5 by default)
  real log_T0_prior;
}

transformed data {
  real t0 = t[1];
  vector[S] log_y = log(y);
  real nu = 5;
}

parameters {
  real log_T0;
  real<lower=0> kill;             // -k_d
  real<lower=0> r;
  real<lower=t0, upper=t_n_max> t_n;
  real<lower=0, upper=1> u;       // position of t_e in [t_n, T[S]] (ignored if held = 0)
  real log_b;
  real<lower=0> sigma;
}

transformed parameters {
  real t_e = held == 1 ? t_n + u * (t[S] - t_n) : t_n;
  vector[S] log_mu;
  for (s in 1:S) {
    real log_T = log_T0 - kill * (fmin(t[s], t_n) - t0) + r * fmax(t[s] - t_e, 0);
    log_mu[s] = log_sum_exp(log_b, log_T);
  }
}

model {
  log_T0 ~ normal(log_T0_prior, 2);
  kill ~ lognormal(log(0.03), 1);        // half-life ~ 23 d, wide
  r ~ lognormal(log(0.01), 1);           // doubling ~ 70 d, wide
  log_b ~ normal(log_b_prior, log_b_sd);
  sigma ~ normal(0, 0.3);
  u ~ uniform(0, 1);
  for (s in 1:S) target += w[s] * student_t_lpdf(log_y[s] | nu, log_mu[s], sigma);
}

generated quantities {
  vector[S] log_lik;
  real held_days = t_e - t_n;
  // depth of the nadir below background (log10), and time below background
  real log10_nadir_vs_b = (log_T0 - kill * (t_n - t0) - log_b) / log(10);
  for (s in 1:S) log_lik[s] = student_t_lpdf(log_y[s] | nu, log_mu[s], sigma);
}
