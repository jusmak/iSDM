Inference_HB <- function(wd, training_data, offset_full) {
  
  #parallelize computation
  options(mc.cores = parallel::detectCores())
  
  # prior for the variation of the linear weights
  linear_sigma = 10
  f_data <- list("y" = training_data$response/training_data$weights, "x" = training_data$covariates,
                 "N" = length(training_data$response), "linear_sigma" = linear_sigma,
                 "weights" = training_data$weights, offset = log(offset_full[,1]) + log(offset_full[,2]))
  
  inits <- list(.1, c(rep(.1,8)))
  names(inits) <- c("alpha", "beta")
  inits_chain1 <- list(inits)
  
  fit_hb_1 <- stan(file=paste(wd, 'R_code/PPP_covariate.stan', sep = '/'), data=f_data, iter=5000, warmup=2000,
                   init = inits_chain1, chains=1, seed=2, refresh=1000, algorithm="NUTS",
                   control=list(adapt_delta=.97, max_treedepth = 15))
  
  return(fit_hb_1)
}