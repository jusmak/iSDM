Inference_glmnet <- function(wd, training_data, offset_full) {
  
  offset <- log(offset_full[,1]) + log(offset_full[,2])
  resp <- training_data$response/training_data$weights
  x_temp <- model.matrix(resp ~ poly(training_data$covariates, degree = 2, raw = TRUE))
  
  fit.glmnet.1 <- glmnet::glmnet(x = x_temp,
                         y = resp,
                         family = "poisson",
                         offset = offset,
                         weights = training_data$weights)

  return(fit.glmnet.1)
}