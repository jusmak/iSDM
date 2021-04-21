data {
  int<lower=1> N_inv;
  int<lower=0> y_inv[N_inv];
  matrix[N_inv,9] x_inv;
  real linear_sigma;
}

parameters {
  vector[9] beta;
}

model {
  // prior models
  to_vector(beta) ~ normal(0,sqrt(linear_sigma));

  // observation model
  for(j in 1:N_inv){
    target += bernoulli_lpmf(y_inv[j] | 1-exp(-exp(x_inv[j,]*beta)));
  }
}
