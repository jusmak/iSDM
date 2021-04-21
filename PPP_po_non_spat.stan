data {
  int<lower=1> N;
  int<lower=0> y[N];
  matrix[N,9] x;
  vector[N] weights_PO;
  real linear_sigma;
}

parameters {
  vector[9] beta;
}

model {
  // prior models
  to_vector(beta) ~ normal(0,sqrt(linear_sigma));

  // observation model
  for(i in 1:N){
    target += poisson_log_lpmf(y[i] | x[i,]*beta + log(weights_PO[i]));
  }
}
