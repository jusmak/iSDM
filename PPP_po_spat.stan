data {
  int<lower=1> N;
  matrix[N,9] x;
  vector[2] coord[N];
  int<lower=0> y[N];
  vector[N] weights_PO;
  real linear_sigma;
}

parameters {
  real<lower=0> lengthscale;
  real<lower=0> sigma;
  vector[9] beta;
  vector[N] eta;
}

transformed parameters {
  // spatial random effect
  vector[N] f;
  matrix[N,N] L_cov;
  matrix[N, N] cov = cov_exp_quad(coord, sigma, lengthscale) + diag_matrix(rep_vector(1e-6, N));
  L_cov = cholesky_decompose(cov);
  f = L_cov * eta;
}

model {
  
  // prior models
  to_vector(beta) ~ normal(0,sqrt(linear_sigma));
  lengthscale ~ normal(0,2);
  sigma ~ normal(0,2);
  eta ~ normal(0,1);
  
  // observation model
  for(i in 1:N){
    target += poisson_log_lpmf(y[i] | x[i,]*beta + f[i] + log(weights_PO[i]));
  }
}
