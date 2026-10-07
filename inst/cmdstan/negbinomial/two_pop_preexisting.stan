functions {
  // Populations on the log scale: exp() of the exponents overflows on long
  // follow-up. An absent population contributes nothing (log 0); the assay
  // background keeps the mean positive, so a small count after clearance does
  // not require a population to explain it.
  real log_ns(real t, real rho_s, real t_end) {
    return t < t_end ? -rho_s * (t - t_end) : negative_infinity();
  }

  real log_nr(real t, real rho_r, real t0_r) {
    return t >= t0_r ? rho_r * (t - t0_r) : negative_infinity();
  }
}

data {
  int<lower=1> S; // Number of steps
  array[S] int<lower=0> N; // observations
  array[S] real T;         // observations
  int<lower=0,upper=1> prior_only;
  real<lower=0> rho_r_prior_sd;   // sd of the half-normal prior on rho_r (per day)
  real<lower=0> rho_s_prior_sd;   // sd of the half-normal prior on rho_s (per day)
  int<lower=0,upper=1> fit_background; // 1: background is sampled; 0: fixed
  real<lower=0> background_fixed;      // background (counts) when not sampled
  real<lower=0> background_prior_mean; // mean of its exponential prior when sampled
}

parameters {
  real<lower=0> rho_r;        // Parameter rho_r (rate for recover)
  real<lower=0> rho_s;        // Parameter rho_s (rate for decay)
  real<upper=T[1]> t0_r;                   // Parameter t_r (time shift)
  real<lower=T[1]> t_end;     // Parameter t_end (time shift)
  real<lower=0> phi;           // NB overdispersion
  array[fit_background] real<lower=0> background_par;
}

transformed parameters {
  // expected count of the assay with no tumour present
  real background = fit_background ? background_par[1] : background_fixed;
}

model {
  vector[S] log_mu;           // log expected values for N

  // Priors
  rho_r ~ normal(0, rho_r_prior_sd);
  rho_s ~ normal(0, rho_s_prior_sd);
  t0_r ~ normal(T[1], T[S] - T[1]);         // Prior for t_r
  t_end ~ normal(T[S], T[S] - T[1]);
  phi ~ gamma(2, 0.1);
  background_par ~ exponential(1 / background_prior_mean);

  for (i in 1:S) {
    log_mu[i] = log_sum_exp({log(background), log_ns(T[i], rho_s, t_end), log_nr(T[i], rho_r, t0_r)});
  }
  if (prior_only == 0) {
    N ~ neg_binomial_2_log(log_mu, phi);
  }
}

generated quantities {
  vector[S] log_lik;             // Log-likelihood for each observation
  vector[S] ns;
  vector[S] nr;
  vector[S] yrep;               // Expected values for N given x

  for (i in 1:S) {
    real log_mu = log_sum_exp({log(background), log_ns(T[i], rho_s, t_end), log_nr(T[i], rho_r, t0_r)});
    ns[i] = exp(log_ns(T[i], rho_s, t_end));
    nr[i] = exp(log_nr(T[i], rho_r, t0_r));
    yrep[i] = exp(log_mu);
    log_lik[i] = neg_binomial_2_log_lpmf(N[i] | log_mu, phi);
  }
}
