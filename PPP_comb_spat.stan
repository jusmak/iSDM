data {
  int<lower=1> N;
  int<lower=1> N_inv;
  matrix[N,9] x;
  matrix[N_inv,9] x_inv;
  vector[2] coord[N];
  int<lower=0> y[N];
  int<lower=0> y_inv[N_inv];
  vector[N] weights_PO;
  real linear_sigma;
}

parameters {
  real<lower=0> lengthscale;
  real<lower=0> sigma;
  real alpha;
  real alpha_bias;
  vector[9] beta;
  vector[N] eta;
}

model {
  // spatial random effect
  vector[N] f;
  {
    matrix[N,N] L_cov;
    matrix[N, N] N_cov = cov_exp_quad(coord, sigma, lengthscale);
    // diagonal elements
    for (n in 1:N) {
      N_cov[n, n] = N_cov[n, n] + 1e-6;
    }
    
    L_cov = cholesky_decompose(N_cov);
    f = L_cov * eta;
  }
  
  // prior models
  to_vector(beta) ~ normal(0,sqrt(linear_sigma));
  alpha_bias ~ normal(0,sqrt(10));
  lengthscale ~ inv_gamma(4,250);
  sigma ~ normal(0,2);
  to_vector(eta) ~ normal(0,1);
  
  // observation model
  for(i in 1:N){
    target += poisson_log_lpmf(y[i] | x[i,]*beta + f[i] + log(weights_PO[i]));
  }
  for(j in 1:N_inv){
    target += bernoulli_lpmf(y_inv[j] | 1-exp(-exp(alpha_bias + x_inv[j,]*beta)));
  }
}
