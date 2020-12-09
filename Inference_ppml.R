Inference_ppml <- function(wd, training_data, offset_full_temp) {
  
  form = ~ poly(cov1, degree = 2, raw = TRUE) + poly(cov2, degree = 2, raw = TRUE) + 
    poly(cov3, degree = 2, raw = TRUE) + poly(cov4, degree = 2, raw = TRUE)
  
  #set data as a data frame in the ppmlasso format
  #coordinates, covariates, response and weights
  #weights are products of the areas and the offsets
  ppm_data <- data.frame(cbind(training_data$coordinates[,1], training_data$coordinates[,2], training_data$covariates[,1],
                   training_data$covariates[,2], training_data$covariates[,3], training_data$covariates[,4],
                   training_data$response, training_data$weights*offset_full_temp[,1]*offset_full_temp[,2]))
  
  colnames(ppm_data) <- c('X', 'Y', colnames(training_data$covariates), 'Pres', 'wt')
  
  fit_ppml = ppmlasso(form, data = ppm_data, n.fits = 50, max.it = 10)
  
  return(fit_ppml)
}