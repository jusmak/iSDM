data {
  int<lower=1> N;
  int<lower=1> N_inv;
  int<lower=0> y[N];
  int<lower=0> y_inv[N_inv];
  matrix[N,9] x;
  matrix[N_inv,9] x_inv;
  vector[N] weights_PO;
  real linear_sigma;
}

parameters {
  real alpha_bias;
  vector[9] beta;
}

model {
  // prior models
  to_vector(beta) ~ normal(0,sqrt(linear_sigma));
  alpha_bias ~ normal(0,sqrt(10));
  
  // observation model
  for(i in 1:N){
    target += poisson_log_lpmf(y[i] | x[i,]*beta + log(weights_PO[i]));
  }
  for(j in 1:N_inv){
    target += bernoulli_lpmf(y_inv[j] | 1-exp(-exp(alpha_bias + x_inv[j,]*beta)));
  }
}
