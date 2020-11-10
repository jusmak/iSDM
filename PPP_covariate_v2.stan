data {
  int<lower=1> N;
  matrix[N,4] x;
  int<lower=0> y[N];
  vector[N] weights;
  int<lower=1> linear_sigma;
  vector[N] expert_p;
}

parameters {
  real alpha;
  vector[4] beta;
}

//transformed parameters {
//  vector[N] log_lambda;
//  log_lambda = alpha + x*beta;
//}

model {
  // prior models
  beta ~ multi_normal(rep_vector(0,4),diag_matrix(rep_vector(sqrt(linear_sigma),4)));
  alpha ~ normal(0,sqrt(10));
  
  // observation model
  //y ~ poisson_log(log_lambda);
  for(i in 1:N){
    //target += poisson_log_lpmf(y[i] | log_lambda[i]) * weights[i];
    target += poisson_log_lpmf(y[i] | alpha + x[i,]*beta) * weights[i] * expert_p[i];
  }
}

//generated quantities {
//  vector[N] log_lik;
//  for (i in 1:N){
//    log_lik[i] = poisson_log_lpmf(y[i] | alpha + x[i,]*beta) * weights[i];
//  }
//}
