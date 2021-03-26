Inference_HB_comb <- function(wd, training_data, PA_data, offset_PO, offset_PA) {
  
  #parallelize computation
  options(mc.cores = parallel::detectCores())
  
  # prior for the variation of the linear weights
  linear_sigma = 10
  #take 
  x_temp <- cbind(training_data$covariates,training_data$covariates^2)
  x_temp_inv <- cbind(PA_data$covariates,PA_data$covariates^2)
  
  f_data <- list("y" = training_data$response/training_data$weights,
                 "y_inv" = PA_data$response, "x" = x_temp, "x_inv" = x_temp_inv,
                 "N" = length(training_data$response), "N_inv" = length(PA_data$response),
                 "linear_sigma" = linear_sigma,
                 "weights_PO" = training_data$weights*training_data$weights_lik,
                 "weights_PA" = PA_data$weights,
                 "offset" = log(offset_PO[,1]) + log(offset_PO[,2]),
                 "offset_inv" = log(offset_PA[,1]) + log(offset_PA[,2]))

  #set initial values
  inits_chain1 <- list(alpha = -1, alpha_bias= -1, beta = rep(-1,8))
  
  #setup the model
  m_hb_1 <- stan_model(file=paste(wd, 'R_code/PPP_covariate_comb_v2.stan', sep = '/'))
  
  #take point estimates of the parameters
  fit_hb_1 <- optimizing(object = m_hb_1, data=f_data, iter = 2000, hessian = TRUE, 
                         seed = 2, init = inits_chain1,
                         draws = 1000, importance_resampling = TRUE, verbose = TRUE)
  
  #check convergenge
  coeff_temp <- fit_hb_1$par
  names(coeff_temp) <- NULL
  conv_ind <- ifelse(all(coeff_temp == -1),F,T)
  
  fit_hb <- list(coefficient = fit_hb_1$par, lik = fit_hb_1$value, convergence_ind = conv_ind)
  return(fit_hb)
}