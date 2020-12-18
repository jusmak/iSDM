Inference_glmnet_ind_validation <-function(wd, model_fits, pred_PA_data) {
  
  
  ind_validation_glmnet_pred_dens <- matrix(NA, ncol = length(pred_PA_data$offset), nrow = length(pred_PA_data$offset))
  ind_validation_glmnet_auc <- matrix(NA, ncol = length(pred_PA_data$offset), nrow = length(pred_PA_data$offset))
  
  for (i in 1:length(pred_PA_data$offset)) {
    for (k in 1:length(pred_PA_data$offset)) {
      
      #define offset for predicting
      offset_temp <- log(pred_PA_data$offset[[k]][,1]) + log(pred_PA_data$offset[[k]][,2])
      
      #define prediction matrix
      x_pred_temp <- model.matrix(~ poly(pred_PA_data$covariates, degree = 2, raw = TRUE))
      
      pred_temp <- predict(object = model_fits$glmnet_model[[i]], newx = x_pred_temp,
                           newoffset = offset_temp)
      pred_lik <- dpois(pred_PA_data$response, exp(pred_temp[,ncol(pred_temp)]), log = TRUE)
      ind_validation_glmnet_pred_dens[i,k] <- sum(pred_lik)
      
      #AUC
      occ_prob <- 1-exp(-exp(pred_temp[,ncol(pred_temp)]))
      suppressMessages(
        myRoc <- pROC::roc(response = pred_PA_data$response, predictor = occ_prob))
      ind_validation_glmnet_auc[i,k] <- myRoc$auc
      
    }
  }
  glmnet_ind_valid <- list(ind_validation_glmnet_pred_dens, ind_validation_glmnet_auc)
  names(glmnet_ind_valid) <- c('ind_validation_glmnet_pred_dens', 'ind_validation_glmnet_auc')
  return(glmnet_ind_valid)
}