Inference_HB_ind_validation <-function(wd, model_fits, pred_PA_data) {
  
  ind_validation_HB_pred_dens <- matrix(NA, ncol = length(model_fits$offset), nrow = length(model_fits$offset))
  ind_validation_HB_auc <- matrix(NA, ncol = length(model_fits$offset), nrow = length(model_fits$offset))
  pred_matrix <- cbind(pred_PA_data$covariates, pred_PA_data$covariates^2)
  
  #loop over different models (offsets used for model fitting)
  for (i in 1:length(model_fits$offset)) {
    beta_samp <- extract(model_fits$HB_model[[i]])
    pred_validation <- matrix(rep(beta_samp$alpha, nrow(pred_matrix)), nrow = nrow(pred_matrix), byrow = TRUE) + pred_matrix %*% t(beta_samp$beta)
    
    #loop over different offset choices for predictions
    for (k in 1:length(model_fits$offset)) {
      pred_lik_pred_dens <- matrix(NA, nrow = nrow(pred_validation), ncol = nrow(extract(model_fits$HB_model[[1]])[[1]]))
      pred_lik_auc <- matrix(NA, ncol = nrow(extract(model_fits$HB_model[[1]])[[1]]))
      for (j in 1:ncol(pred_validation)) {
        pred_lik_pred_dens[,j] <- dpois(pred_PA_data$response, exp(pred_PA_data$offset[[k]] + pred_validation[,j]), log = TRUE)
        
        occ_prob <- 1-exp(-exp(pred_PA_data$offset[[k]] + pred_validation[,j]))
        suppressMessages(
          myRoc <- pROC::roc(response = pred_PA_data$response, predictor = occ_prob))          
        pred_lik_auc[j] <- myRoc$auc
      }
      pred_dens <- apply(pred_lik_pred_dens, 1, function(x) -log(ncol(pred_lik_pred_dens)) + log(sum(exp(x))))
      ind_validation_HB_pred_dens[i,k] <- sum(pred_dens)
      ind_validation_HB_auc[i,k] <- mean(pred_lik_auc)
    }
  }
  HB_ind_valid <- list(ind_validation_HB_pred_dens, ind_validation_HB_auc)
  names(HB_ind_valid) <- c('ind_validation_HB_pred_dens', 'ind_validation_HB_auc')
  return(HB_ind_valid)
}

