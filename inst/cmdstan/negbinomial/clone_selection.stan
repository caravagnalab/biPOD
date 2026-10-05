// Joint piecewise-exponential model of clone-specific burden.
//
// biPOD-style exponential growth with one breakpoint t_b shared by all clones
// of a series (response -> regrowth). Each clone k has its own
// response-phase rate rho1[k] and regrowth-phase rate rho2[k]:
//
//   log mu[k, s] = a[k] + rho1[k] * (min(T[s], t_b) - T[1])
//                       + rho2[k] * max(T[s] - t_b, 0) + depth_offset[k, s]
//
// Observations are alt reads over the clone-specific SVs, with a depth offset
// and negative-binomial noise. The initial size a[k] is free, so clones that
// are undetected at baseline are allowed.

data {
  int<lower=1> K;                     // clones
  int<lower=2> S;                     // timepoints
  array[S] real T;                    // sorted sampling times (days)
  array[K, S] int<lower=0> y;         // alt reads per clone and timepoint
  matrix[K, S] depth_offset;          // log(depth / reference depth)
  vector[K] a_prior;                  // prior centre for log initial burden
  real t_b_prior;                     // prior centre for the shared breakpoint
  real<lower=0> t_b_prior_sd;
  real<lower=0> rho_prior_sd;
}

parameters {
  vector[K] a;
  vector[K] rho1;
  vector[K] rho2;
  real<lower=T[1], upper=T[S]> t_b;
  real<lower=0> phi;
}

transformed parameters {
  matrix[K, S] log_mu;
  for (k in 1:K) {
    for (s in 1:S) {
      log_mu[k, s] = a[k] + rho1[k] * (fmin(T[s], t_b) - T[1])
                     + rho2[k] * fmax(T[s] - t_b, 0) + depth_offset[k, s];
    }
  }
}

model {
  a ~ normal(a_prior, 3);
  rho1 ~ normal(0, rho_prior_sd);
  rho2 ~ normal(0, rho_prior_sd);
  t_b ~ normal(t_b_prior, t_b_prior_sd);
  phi ~ gamma(2, 0.1);

  for (k in 1:K) y[k] ~ neg_binomial_2_log(log_mu[k], phi);
}

generated quantities {
  // yrep is -1 where the posterior draw implies an overflowing mean
  array[K, S] int yrep;
  matrix[K, S] log_lik;
  // log fold change of each clone's burden over the whole window
  vector[K] log_fc;
  // change in log relative frequency over the window (> 0: positive selection)
  vector[K] delta_log_freq;

  {
    vector[K] log_start = a;
    vector[K] log_end;
    for (k in 1:K) {
      for (s in 1:S) {
        yrep[k, s] = log_mu[k, s] < 20 ? neg_binomial_2_log_rng(log_mu[k, s], phi) : -1;
        log_lik[k, s] = neg_binomial_2_log_lpmf(y[k, s] | log_mu[k, s], phi);
      }
      log_fc[k] = rho1[k] * (t_b - T[1]) + rho2[k] * (T[S] - t_b);
      log_end[k] = a[k] + log_fc[k];
    }
    delta_log_freq = (log_end - log_sum_exp(log_end)) - (log_start - log_sum_exp(log_start));
  }
}
