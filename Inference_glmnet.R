Inference_glmnet <- function(wd, training_data, offset_full) {
  
  form = ~ poly(cov1, degree = 2, raw = TRUE) + poly(cov2, degree = 2, raw = TRUE) + 
    poly(cov3, degree = 2, raw = TRUE) + poly(cov4, degree = 2, raw = TRUE)
  fit.glmnet.1 <- glmnet(formula = form, x = training_data$covariates, 
                         y = training_data$response, family = "poisson", weights = apply(offset_full,1,prod))

  return(fit.glmnet.1)
}