Inference_glm <- function(wd, training_data, offset_temp) {
  
  offset_temp <- log(offset_temp[,1]) - mean(log(offset_temp[,1])) + log(offset_temp[,2]) - mean(log(offset_temp[,2]))
  form <- y ~ cov1 + I(cov1^2) + cov2 + I(cov2^2) + cov3 + I(cov3^2) + cov4 + I(cov4^2)
  data_temp <- cbind.data.frame(training_data$response/training_data$weights, training_data$covariates)
  names(data_temp) <- c('y', 'cov1', 'cov2', 'cov3', 'cov4')
  fit.glm.1 <- glm(formula = form, data = data_temp, family = "poisson",
                   offset = offset_temp,
                   weights = training_data$weights,
                   maxit=100)

  return(fit.glm.1)
}