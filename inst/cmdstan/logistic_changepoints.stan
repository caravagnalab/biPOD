
functions {

  real mean_t (int n0, real t0, real t, real K, real[] changing_times, real[] rho_array) {
    real nj = n0;
    real res = n0;
    int J = num_elements(changing_times);

    if (t <= changing_times[1]) {
      return nj * K / (nj + (K - nj) * exp(-rho_array[1] * (t - t0)));
    } else {
      nj = nj * K / (nj + (K - nj) * exp(-rho_array[1] * (t - t0)));

      for (j in 2:J) {
        if (t >= changing_times[j]) {
          nj = nj * K / (nj + (K - nj) * exp(-rho_array[j] * (changing_times[j] - changing_times[j-1])));
        } else {
          res = nj * K / (nj + (K - nj) * exp(-rho_array[j] * (t - changing_times[j-1])));
          return(res);
        }
      }
      res = nj * K / (nj + (K - nj) * exp(-rho_array[J+1] * (t - changing_times[J])));
      return(res);
    }
  }
}

data {
  int<lower=1> S; // Number of steps
  int<lower=1> G; // Number of wondows

  int <lower=0> N[S];      // observations
  real T[S];      // observations

  real<lower=0> prior_K;
  real changing_times_prior[G - 1];
  real dt;
}

parameters {
  real rho[G];
  real<lower=-dt, upper=dt> changing_times_unit[G-1];
  real<lower=prior_K> K; // carrying capacity
}

transformed parameters {
  real changing_times[G-1];

  for (i in 2:G) {
    changing_times[i-1] = changing_times_prior[i-1] + changing_times_unit[i-1];
  }
}

model {
  target += normal_lpdf(K | prior_K, prior_K); // sample the carrying capacity

  for (i in 1:G) {
    target += normal_lpdf(rho[i] | 0, 1);
  }

  for (i in 2:G) {
    target += normal_lpdf(changing_times_unit[i-1] | 0, dt);
  }

  for (i in 1:S) {
    target += poisson_lpmf(N[i] | mean_t(N[1], T[1], T[i], K, changing_times, rho));
  }
}
