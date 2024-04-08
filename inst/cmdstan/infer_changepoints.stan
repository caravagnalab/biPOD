
functions {

  real mean_t (int n0, real t0, real t, vector changing_times, vector rho_array) {
    real res = n0;
    int J = num_elements(changing_times);

    if (t <= changing_times[1]) {
      return n0 * exp(rho_array[1] * (t - t0));
    } else {
      res = res * exp(rho_array[1] * (changing_times[1] - t0));
      for (j in 2:J) {
        if (t >= changing_times[j]) {
          res = res * exp(rho_array[j] * (changing_times[j] - changing_times[j-1]));
        } else {
          res = res * exp(rho_array[j] * (t - changing_times[j-1]));
          return(res);
        }
      }
      res = res * exp(rho_array[J+1] * (t - changing_times[J]));
      return(res);
    }
  }
}

data {
  int<lower=1> S; // Number of steps
  int<lower=1> G; // Number of wondows

  array[S] int<lower=0> N; // observations
  array[S] real T;         // observations

  array[G-1] real changing_times_prior;
  real dt;
}

parameters {
  vector[G] rho;
  array[G-1] real<lower=-dt, upper=dt> changing_times_unit;
}

transformed parameters {
  vector[G-1] changing_times;

  for (i in 2:G) {
    changing_times[i-1] = changing_times_prior[i-1] + changing_times_unit[i-1];
  }
}

model {
  for (i in 1:G) {
    target += normal_lpdf(rho[i] | 0, 1);
  }

  for (i in 2:G) {
    target += normal_lpdf(changing_times_unit[i-1] | 0, dt);
  }

  for (i in 1:S) {
    target += poisson_lpmf(N[i] | mean_t(N[1], T[1], T[i], changing_times, rho));
  }
}
