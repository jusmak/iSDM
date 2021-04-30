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


model {
  
  // spatial random effect
  vector[N] f;
  vector[N] f_nu;
  {
    matrix[N,N] L_cov;
    matrix[N, N] N_cov = cov_exp_quad(coord, sigma, lengthscale);
    // diagonal elements
    for (n in 1:N) {
      N_cov[n, n] = N_cov[n, n] + 1e-6;
    }
    
    L_cov = cholesky_decompose(N_cov);
    f = L_cov * eta;
    f_nu = f - L_cov * x / (x' * L_cov * x) * x' * f;
  }
  
  // prior models
  to_vector(beta) ~ normal(0,sqrt(linear_sigma));
  lengthscale ~ inv_gamma(4,250);
  sigma ~ normal(0,2);
  to_vector(eta) ~ normal(0,1);
  
  // observation model
  for(i in 1:N){
    target += poisson_log_lpmf(y[i] | x[i,]*beta + f[i] + log(weights_PO[i]));
  }
}
