functions {
  real integrated_r(real t, real t0, vector t_array, vector rho_array) {
    int n_t = num_elements(t_array);
    real res = 0;

    if (n_t == 0) return rho_array[1] * (t - t0);
    if (t <= t_array[1]) return rho_array[1] * (t - t0);

    res = rho_array[1] * (t_array[1] - t0);
    for (i in 2:n_t) {
      if (t <= t_array[i])
        return res + rho_array[i] * (t - t_array[i - 1]);
      res += rho_array[i] * (t_array[i] - t_array[i - 1]);
    }
    return res + rho_array[n_t + 1] * (t - t_array[n_t]);
  }

  real mean_t(real t, real t0, vector t_array, vector rho_array, real K) {
    // y(t) = K - (K - 1) * exp(-∫ rho)
    return K - (K - 1) * exp(-integrated_r(t, t0, t_array, rho_array));
  }
}

data {
  int<lower=1> S;             // # observations
  int<lower=1> G;             // # segments (=> length(rho) = G)
  array[S] int<lower=0> N;    // counts
  array[S] real T;            // times
  vector[G - 1] t_array;      // breakpoints (strictly increasing)
  int<lower=0,upper=1> prior_only;
}

parameters {
  vector[G] rho;              // segment rates
  real<upper=T[1]> t0;        // initiation time
  real<lower=0> sigma;        // log-normal noise sd
  real<lower=1> K;            // asymptote
}

model {
  // Priors (adjust to taste)
  rho ~ normal(0, 1);
  t0  ~ normal(T[1], 30);
  sigma ~ exponential(1);
  K ~ normal(max(N), max(N));

  if (prior_only == 0) {
    vector[S] mu_log; // location for LogNormal (mean on log scale)
    for (i in 1:S) {
      real mu_pred = mean_t(T[i], t0, t_array, rho, K);
      mu_log[i] = log(mu_pred);
    }
    N ~ lognormal(mu_log, sigma);
  }
}

generated quantities {
  vector[S] log_lik;
  array[S] real yrep;
  vector[S] mu_pred;

  if (prior_only == 0) {
    for (i in 1:S) {
      mu_pred[i] = mean_t(T[i], t0, t_array, rho, K);
      log_lik[i] = lognormal_lpdf(N[i] | log(mu_pred[i]), sigma);
      yrep[i]    = lognormal_rng(log(mu_pred[i]), sigma);
    }
  }
}
