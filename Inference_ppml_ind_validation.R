Inference_ppml_ind_validation <-function(wd, model_fits, pred_PA_data) {
  
  
  ind_validation_ppml_pred_dens <- matrix(NA, ncol = length(pred_PA_data$offset), nrow = length(pred_PA_data$offset))
  ind_validation_ppml_auc <- matrix(NA, ncol = length(pred_PA_data$offset), nrow = length(pred_PA_data$offset))
  
  for (i in 1:length(pred_PA_data$offset)) {
    for (k in 1:length(pred_PA_data$offset)) {
      
      #define data for predicting
      ppm_data_pred <- data.frame(pred_PA_data$coordinates, pred_PA_data$covariates,
                                  pred_PA_data$offset[[k]][,1] * pred_PA_data$offset[[k]][,2])
      
      colnames(ppm_data_pred) <- c('X', 'Y', colnames(training_data$covariates), 'wt')
      
      pred_temp <- predict(model_fits$ppml_model[[i]], newdata = ppm_data_pred)
      pred_lik <- dpois(pred_PA_data$response, exp(pred_temp[,ncol(pred_temp)]), log = TRUE)
      ind_validation_ppml_pred_dens[i,k] <- sum(pred_lik)
      
      #AUC
      occ_prob <- 1-exp(-exp(pred_temp[,ncol(pred_temp)]))
      suppressMessages(
        myRoc <- pROC::roc(response = pred_PA_data$response, predictor = occ_prob))
      ind_validation_ppml_auc[i,k] <- myRoc$auc
      
    }
  }
  ppml_ind_valid <- list(ind_validation_ppml_pred_dens, ind_validation_ppml_auc)
  names(ppml_ind_valid) <- c('ind_validation_ppml_pred_dens', 'ind_validation_ppml_auc')
  return(ppml_ind_valid)
}