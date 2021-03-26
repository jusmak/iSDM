Inference_HB_comb <- function(wd, training_data, PA_data, offset_PO, offset_PA, weights_temp_set) {
  
  #parallelize computation
  options(mc.cores = parallel::detectCores())
  
  # prior for the variation of the linear weights
  linear_sigma = 10
  #take 
  x_temp <- cbind(training_data$covariates,training_data$covariates^2)
  x_temp_inv <- cbind(PA_data$covariates,PA_data$covariates^2)
  
  #move presence and their related background points locations slightly
  training_data$coordinates[training_data$response==1,1] <- training_data$coordinates[training_data$response==1,1] -.25
  training_data$coordinates[training_data$response==0 & training_data$weights==.5,1] <- training_data$coordinates[training_data$response==0 & training_data$weights==.5,1] +.25

  f_data <- list("y" = training_data$response/training_data$weights,
                 "y_inv" = PA_data$response, "x" = x_temp, "x_inv" = x_temp_inv,
                 "coord" = training_data$coordinates,
                 "N" = length(training_data$response), "N_inv" = length(PA_data$response),
                 "linear_sigma" = linear_sigma,
                 "weights_PO" = training_data$weights*training_data$weights_lik,
                 "weights_PA" = PA_data$weights,
                 "offset" = log(offset_PO[,1]) + log(offset_PO[,2]),
                 "offset_inv" = log(offset_PA[,1]) + log(offset_PA[,2]))

  
  
  if (weights_temp_set==1) {
    #set initial values
    inits_chain1 <- list(alpha = -1, alpha_bias= -1, beta = rep(-1,8))
    
    #setup the model
    m_hb_1 <- paste(wd, 'R_code/PPP_covariate_comb_v2.stan', sep = '/')
  } else {
    #set initial values
    inits_chain1 <- list(lengthscale = .1, sigma = .1, alpha = -1, alpha_bias= -1, beta = rep(-1,8), eta = rep(.1,length(training_data$response)))
    
    #setup the model
    m_hb_1 <- paste(wd, 'R_code/PPP_covariate_comb_v3.stan', sep = '/')  
  }
  
  
  #take point estimates of the parameters
#  fit_hb_1 <- optimizing(object = m_hb_1, data=f_data, iter = 2000, hessian = TRUE, 
#                         seed = 2, init = inits_chain1,
#                         draws = 1000, importance_resampling = TRUE, verbose = TRUE)
  fit_hb_1 <- stan(file=m_hb_1, data=f_data, iter=5000, warmup=2000,
                   init = inits_chain1, chains=1, seed=2, refresh=1000, algorithm="NUTS",
                   control=list(adapt_delta=.97, max_treedepth = 15))
  
  fit_hb <- fit_hb_1
  #check convergenge
#  if (weights_temp_set == 1) {
#    conv_ind <- ifelse(all.equal(fit_hb_1$par,-.1),F,T)
#    fit_hb <- list(coefficient = fit_hb_1$par, lik = fit_hb_1$value, convergence_ind = conv_ind)
#    
#  } else {
#    temp_parameters <- fit_hb_1$par[3:12]
#    conv_ind <- ifelse(all.equal(temp_parameters,-.1),F,T)
#    store_parameters <- fit_hb_1$par[1:2]
#    fit_hb <- list(coefficient = temp_parameters, gp_parameters = store_parameters, lik = fit_hb_1$value, convergence_ind = conv_ind)
#  }
  
  return(fit_hb)
}