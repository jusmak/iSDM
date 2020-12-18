Inference_ppml_interp <-function(wd, model_fits, training_data, offset_full) {
  
  
  train_validation_glmnet_pred_dens <- matrix(NA, ncol = length(offset_full), nrow = length(offset_full))
  train_validation_glmnet_auc <- matrix(NA, ncol = length(offset_full), nrow = length(offset_full))
  
  for (i in 1) {
    for (k in 1) {
      
      #define offset for predicting
      offset_temp <- log(offset_full[[k]][,1]) + log(offset_full[[k]][,2])
      
      #define prediction matrix
      x_pred_temp <- model.matrix(~ poly(training_data$covariates, degree = 2, raw = TRUE))
      
      pred_temp <- predict(object = model_fits$glmnet_model[[i]], newx = x_pred_temp,
                           newoffset = offset_temp)
      pred_lik <- dpois(training_data$response, exp(pred_temp[,ncol(pred_temp)]), log = TRUE)
      train_validation_glmnet_pred_dens[i,k] <- sum(pred_lik)
      
      #AUC
      occ_prob <- 1-exp(-exp(pred_temp[,ncol(pred_temp)]))
      suppressMessages(
        myRoc <- pROC::roc(response = training_data$response, predictor = occ_prob))
      train_validation_glmnet_auc[i,k] <- myRoc$auc
      
    }
  }
  glmnet_train_valid <- list(train_validation_glmnet_pred_dens, train_validation_glmnet_auc)
  names(glmnet_train_valid) <- c('train_validation_glmnet_pred_dens', 'train_validation_glmnet_auc')
  return(glmnet_train_valid)
}