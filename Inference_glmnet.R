Inference_glmnet <- function(wd, training_data, offset_full) {
  
  offset = log(offset_full[,1]) + log(offset_full[,2])
  form = ~ poly(cov1, degree = 2, raw = TRUE) + poly(cov2, degree = 2, raw = TRUE) + 
    poly(cov3, degree = 2, raw = TRUE) + poly(cov4, degree = 2, raw = TRUE)
  fit.glmnet.1 <- glmnet(formula = form, x = training_data$covariates, 
                         y = training_data$response/training_data$weights, family = "poisson",
                         offset = offset,
                         weights = training_data$weights)

  return(fit.glmnet.1)
}