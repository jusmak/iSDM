data {
  int<lower=1> N;
  matrix[N,4] x;
  int<lower=0> y[N];
  vector[N] weights;
  vector[N] offset;
  int<lower=1> linear_sigma;
}

transformed data {
  matrix[N,4] x_squared = x .* x;
  matrix[N,8] x_full = append_col(x,x_squared);
}

parameters {
  real alpha;
  vector[8] beta;
}

//transformed parameters {
//  vector[N] log_lambda;
//  log_lambda = alpha + x*beta;
//}

model {
  // prior models
  beta ~ multi_normal(rep_vector(0,8),diag_matrix(rep_vector(sqrt(linear_sigma),8)));
  alpha ~ normal(0,sqrt(10));
  
  // observation model
  for(i in 1:N){
    target += poisson_log_lpmf(y[i] | offset[i] + alpha + x_full[i,]*beta) * weights[i];
  }
}

//generated quantities {
//  vector[N] log_lik;
//  for (i in 1:N){
//    log_lik[i] = poisson_lpmf(y[i] | offset[i] + alpha + x[i,]*beta) * weights[i];
//  }
//}
