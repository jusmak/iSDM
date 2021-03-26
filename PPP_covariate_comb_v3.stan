data {
  int<lower=1> N;
  int<lower=1> N_inv;
  matrix[N,8] x;
  matrix[N_inv,8] x_inv;
  vector[2] coord[N];
  int<lower=0> y[N];
  int<lower=0> y_inv[N_inv];
  vector[N] weights_PO;
  vector[N_inv] weights_PA;
  vector[N] offset;
  vector[N_inv] offset_inv;
  real linear_sigma;
}

parameters {
  real<lower=0> lengthscale;
  real<lower=0> sigma;
  real alpha;
  real alpha_bias;
  vector[8] beta;
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
  alpha ~ normal(0,sqrt(10));
  alpha_bias ~ normal(0,sqrt(10));
  lengthscale ~ normal(0,2);
  sigma ~ normal(0,2);
  eta ~ normal(0,1);
  
  // observation model
  for(i in 1:N){
    target += poisson_log_lpmf(y[i] | offset[i] + alpha + x[i,]*beta + f[i]) * weights_PO[i];
  }
  for(j in 1:N_inv){
    target += bernoulli_lpmf(y_inv[j] | 1-exp(-exp(offset_inv[j] + alpha + alpha_bias + x_inv[j,]*beta))) * weights_PA[j];
  }
  
}

//generated quantities {
//  vector[N] log_lik;
//  for (i in 1:N){
//    log_lik[i] = poisson_lpmf(y[i] | offset[i] + alpha + x[i,]*beta) * weights[i];
//  }
//}
