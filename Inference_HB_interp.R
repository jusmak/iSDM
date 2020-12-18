Inference_HB_interp <- function(wd, model_fits, training_data, offset_full, samp_temp) {
  
  #make a prediction in training locations
  train_validation_HB_pred_dens <- matrix(NA, ncol = length(offset_full), nrow = length(offset_full))
  train_validation_HB_auc <- matrix(NA, ncol = length(offset_full), nrow = length(offset_full))
  pred_matrix <- cbind(training_data$covariates, training_data$covariates^2)
  
  #define matrices for storing values

  for (i in 1:length(offset_full)) {
    beta_samp <- extract(model_fits$HB_model[[i]])
    pred_validation <- matrix(rep(beta_samp$alpha, nrow(pred_matrix)), nrow = nrow(pred_matrix), byrow = TRUE) + pred_matrix %*% t(beta_samp$beta)
    
    for (k in 1:length(offset_full)) {

      #define offset for predicting
      offset_temp <- log(offset_full[[k]][,1]) + log(offset_full[[k]][,2])
      
      #predictive density
      pred_lik_pred_dens <- apply(pred_validation[,samp_temp], 2, function(x) dpois(training_data$response,
                                                                        exp(offset_temp + x), log = TRUE))
      #AUC
      get_auc <- function(x) {
        suppressMessages(
        myRoc <- pROC::roc(response = training_data$response, predictor = 1-exp(-exp(offset_temp + x))))
        return(myRoc$auc)
      }
      pred_lik_auc <- apply(pred_validation[,samp_temp], 2, get_auc)
 
      #sum of predictive densities over the training locations
      pred_dens <- apply(pred_lik_pred_dens, 1, function(x) -log(ncol(pred_lik_pred_dens)) + log(sum(exp(x))))
      train_validation_HB_pred_dens[i,k] <- sum(pred_dens)
      #mean AUC over the training locations
      train_validation_HB_auc[i,k] <- mean(pred_lik_auc)
    }
  }
  HB_train_valid <- list(train_validation_HB_pred_dens, train_validation_HB_auc)
  names(HB_train_valid) <- c('train_validation_HB_pred_dens', 'train_validation_HB_auc')
  return(HB_train_valid)
}
