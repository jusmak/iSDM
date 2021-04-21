Inference_non_spat_po <- function(wd, data, n_samples, burn_in) {
  
  #parallelize computation
  options(mc.cores = parallel::detectCores())
  
  # prior for the variation of the linear weights
  linear_sigma = 10
  
  # add a constant parameter in covariates
  x_temp <- cbind(rep(1,length(data$training_data$response)),data$training_data$covariates)
  x_inv_temp <- cbind(rep(1,length(data$inventory_data$response)),data$inventory_data$covariates)
  
  f_data <- list("N" = length(data$training_data$response),
                 "y" = data$training_data$response/data$training_data$weights,
                 "x" = x_temp, "weights_PO" = data$training_data$weights,
                 "linear_sigma" = linear_sigma)
  
  #set initial values
  inits_chain1 <- list(beta = rep(-1,9))
  
  #setup the model
  m_hb_1 <- paste(wd, 'R_code/PPP_po_non_spat.stan', sep = '/')
  
  fit_hb <- stan(file=m_hb_1, data=f_data, iter=n_samples, warmup=burn_in,
                 init = inits_chain1, chains=1, seed=2, refresh=1000, algorithm="NUTS",
                 control=list(adapt_delta=.97, max_treedepth = 15))
  return(fit_hb)
}

Inference_non_spat_pa <- function(wd, data, n_samples, burn_in) {
  
  #parallelize computation
  options(mc.cores = parallel::detectCores())
  
  # prior for the variation of the linear weights
  linear_sigma = 10
  
  # add a constant parameter in covariates
  x_temp <- cbind(rep(1,length(data$training_data$response)),data$training_data$covariates)
  x_inv_temp <- cbind(rep(1,length(data$inventory_data$response)),data$inventory_data$covariates)
  
  f_data <- list("N_inv" = length(data$inventory_data$response),
                 "y_inv" = data$inventory_data$response,
                 "x_inv" = x_inv_temp, "linear_sigma" = linear_sigma)
  
  #set initial values
  inits_chain1 <- list(beta = rep(-1,9))
  
  #setup the model
  m_hb_1 <- paste(wd, 'R_code/PPP_pa_non_spat.stan', sep = '/')
  
  fit_hb <- stan(file=m_hb_1, data=f_data, iter=n_samples, warmup=burn_in,
                 init = inits_chain1, chains=1, seed=2, refresh=1000, algorithm="NUTS",
                 control=list(adapt_delta=.97, max_treedepth = 15))
  return(fit_hb)
}

Inference_spat_po <- function(wd, data, n_samples, burn_in) {
  
  #parallelize computation
  options(mc.cores = parallel::detectCores())
  
  # prior for the variation of the linear weights
  linear_sigma = 10
  
  # add a constant parameter in covariates
  x_temp <- cbind(rep(1,length(data$training_data$response)),data$training_data$covariates)

  f_data <- list("N" = length(data$training_data$response),
                 "y" = data$training_data$response/data$training_data$weights,
                 "x" = x_temp,  "coord" = data$training_data$coordinates,
                 "weights_PO" = data$training_data$weights,
                 "linear_sigma" = linear_sigma)
  
  #set initial values
  inits_chain1 <- list(lengthscale = .1, sigma = .1, beta = rep(-1,9),
                       eta = rep(.1,length(training_data$response)))
  
  #setup the model
  m_hb_1 <- paste(wd, 'R_code/PPP_po_spat.stan', sep = '/')
  
  fit_hb <- stan(file=m_hb_1, data=f_data, iter=n_samples, warmup=burn_in,
                 init = inits_chain1, chains=1, seed=2, refresh=1000, algorithm="NUTS",
                 control=list(adapt_delta=.97, max_treedepth = 15))
  return(fit_hb)
}

Inference_restr_spat_po <- function(wd, data, n_samples, burn_in) {
  
  #parallelize computation
  options(mc.cores = parallel::detectCores())
  
  # prior for the variation of the linear weights
  linear_sigma = 10
  
  # add a constant parameter in covariates
  x_temp <- cbind(rep(1,length(data$training_data$response)),data$training_data$covariates)

  f_data <- list("N" = length(data$training_data$response),
                 "y" = data$training_data$response/data$training_data$weights,
                 "x" = x_temp,  "coord" = data$training_data$coordinates,
                 "weights_PO" = data$training_data$weights,
                 "linear_sigma" = linear_sigma)
  
  #set initial values
  inits_chain1 <- list(lengthscale = .1, sigma = .1, beta = rep(-1,9),
                       eta = rep(.1,length(training_data$response)))
  
  #setup the model
  m_hb_1 <- paste(wd, 'R_code/PPP_po_restr_spat.stan', sep = '/')
  
  fit_hb <- stan(file=m_hb_1, data=f_data, iter=n_samples, warmup=burn_in,
                 init = inits_chain1, chains=1, seed=2, refresh=1000, algorithm="NUTS",
                 control=list(adapt_delta=.97, max_treedepth = 15))
  return(fit_hb)
}


Inference_non_spat_comb <- function(wd, data, n_samples, burn_in) {
  
  #parallelize computation
  options(mc.cores = parallel::detectCores())
  
  # prior for the variation of the linear weights
  linear_sigma = 10
  
  # add a constant parameter in covariates
  x_temp <- cbind(rep(1,length(data$training_data$response)),data$training_data$covariates)
  x_inv_temp <- cbind(rep(1,length(data$inventory_data$response)),data$inventory_data$covariates)
  
  f_data <- list("N" = length(data$training_data$response),
                 "N_inv" = length(data$inventory_data$response),
                 "y" = data$training_data$response/data$training_data$weights,
                 "y_inv" = data$inventory_data$response, "x" = x_temp,
                 "x_inv" = x_inv_temp, "weights_PO" = data$training_data$weights,
                 "linear_sigma" = linear_sigma)
  
  #set initial values
  inits_chain1 <- list(alpha_bias= -1, beta = rep(-1,9))
  
  #setup the model
  m_hb_1 <- paste(wd, 'R_code/PPP_comb_non_spat.stan', sep = '/')
  
  fit_hb <- stan(file=m_hb_1, data=f_data, iter=n_samples, warmup=burn_in,
                 init = inits_chain1, chains=1, seed=2, refresh=1000, algorithm="NUTS",
                 control=list(adapt_delta=.97, max_treedepth = 15))
  return(fit_hb)
}



Inference_spat_comb <- function(wd, data, n_samples, burn_in) {
  
  #parallelize computation
  options(mc.cores = parallel::detectCores())
  
  # prior for the variation of the linear weights
  linear_sigma = 10
  
  # add a constant parameter in covariates
  x_temp <- cbind(rep(1,length(data$training_data$response)),data$training_data$covariates)
  x_inv_temp <- cbind(rep(1,length(data$inventory_data$response)),data$inventory_data$covariates)
  
  f_data <- list("N" = length(data$training_data$response), "N_inv" = length(data$inventory_data$response),
                 "y" = data$training_data$response/data$training_data$weights,
                 "y_inv" = data$inventory_data$response, "x" = x_temp,
                 "x_inv" = x_inv_temp, "coord" = data$training_data$coordinates,
                 "weights_PO" = data$training_data$weights, "linear_sigma" = linear_sigma)
  
  #set initial values
  inits_chain1 <- list(lengthscale = .1, sigma = .1, alpha_bias= -1,
                       beta = rep(-1,9), eta = rep(.1,length(training_data$response)))
  
  #setup the model
  m_hb_1 <- paste(wd, 'R_code/PPP_comb_spat.stan', sep = '/')
  fit_hb <- stan(file=m_hb_1, data=f_data, iter=n_samples, warmup=burn_in,
                 init = inits_chain1, chains=1, seed=2, refresh=1000, algorithm="NUTS",
                 control=list(adapt_delta=.97, max_treedepth = 15))
  return(fit_hb)
}



Inference_restr_spat_comb <- function(wd, data, n_samples, burn_in) {
  
  #parallelize computation
  options(mc.cores = parallel::detectCores())
  
  # prior for the variation of the linear weights
  linear_sigma = 10
  
  # add a constant parameter in covariates
  x_temp <- cbind(rep(1,length(data$training_data$response)),data$training_data$covariates)
  x_inv_temp <- cbind(rep(1,length(data$inventory_data$response)),data$inventory_data$covariates)
  
  f_data <- list("N" = length(data$training_data$response), "N_inv" = length(data$inventory_data$response),
                 "y" = data$training_data$response/data$training_data$weights,
                 "y_inv" = data$inventory_data$response, "x" = x_temp,
                 "x_inv" = x_inv_temp, "coord" = data$training_data$coordinates,
                 "weights_PO" = data$training_data$weights, "linear_sigma" = linear_sigma)
  
  #set initial values
  inits_chain1 <- list(lengthscale = .1, sigma = .1, alpha_bias= -1,
                       beta = rep(-1,9), eta = rep(.1,length(training_data$response)))
  
  #setup the model
  m_hb_1 <- paste(wd, 'R_code/PPP_comb_restr_spat.stan', sep = '/')
  
  fit_hb <- stan(file=m_hb_1, data=f_data, iter=n_samples, warmup=burn_in,
                 init = inits_chain1, chains=1, seed=2, refresh=1000, algorithm="NUTS",
                 control=list(adapt_delta=.97, max_treedepth = 15))
  return(fit_hb)
}

