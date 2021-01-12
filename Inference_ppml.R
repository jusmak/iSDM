Inference_ppml <- function(wd, training_data, offset_full) {
  
  form = ~ cov1 + cov2 + cov3 + cov4 + I(cov1^2) + I(cov2^2) + I(cov3^2) + 
    I(cov4^2) + I(cov1*cov2) + I(cov1*cov3) + I(cov1*cov4) + I(cov2*cov3) +
    I(cov2*cov4) + I(cov3*cov4)
  
  #set data as a data frame in the ppmlasso format
  #coordinates, covariates, response and weights
  ppm_data <- data.frame(cbind(training_data$coordinates, training_data$covariates,
                   training_data$response, training_data$weights))
  
  colnames(ppm_data) <- c('X', 'Y', colnames(training_data$covariates), 'Pres', 'wt')
  
  fit_ppml = ppmlasso(form, data = ppm_data, n.fits = 50, max.it = 10)
  
  fit_temp <- list(form, fit_ppml$beta, fit_ppml$likelihood, fit_ppml$criterion.matrix[,2],
                        fit_ppml$lambdas, fit_ppml$s.means, fit_ppml$s.sds, fit_ppml$lambdas)
  names(fit_temp) <- c('formula', 'coefficients', 'pen_lik', 'BIC', 'lambda', 's_mean', 's_sd')
  
  return(fit_temp)
}