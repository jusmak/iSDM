Plot_thresholds <- function(model_fits, ind_samples){
  
  TSS_thres <- matrix(NA, nrow = length(model_fits$offset))
  TP05_thres <- matrix(NA, nrow = length(model_fits$offset))
  AUC_int <- matrix(NA, nrow = length(model_fits$offset))
  pred_dens_int <- matrix(NA, nrow = length(model_fits$offset))
  
  #make a prediction in training locations
  pred_matrix <- cbind(model_fits$training_data$covariates, model_fits$training_data$covariates^2)
  
  #add offsets
  for (i in 1:length(model_fits$offset)) {
    beta_samp <- extract(model_fits$HB_model[[i]])
    pred <- matrix(rep(beta_samp$alpha[ind_samples], nrow(pred_matrix)), nrow = nrow(pred_matrix), byrow = TRUE) + pred_matrix %*% t(beta_samp$beta[ind_samples,])
    
    offset_temp <- log(model_fits$offset[[i]][,1]) - mean(log(model_fits$offset[[i]][,1])) +
      log(model_fits$offset[[i]][,2]) - mean(log(model_fits$offset[[i]][,2]))
    pred <- exp(pred + offset_temp + log(model_fits$training_data$weights))
    pred_mean <- apply(pred, 1, mean)
    occ_prob_mean <- 1-exp(-(pred_mean))
    roc_temp <- pROC::roc(response = model_fits$training_data$response, predictor = occ_prob_mean)
    TSS_thres[i] <- log(pred)[which.max(roc_temp$sensitivities+roc_temp$specificities-1)]
    TP05_thres[i] <- log(pred)[which(roc_temp$sensitivities<=.95)[1]]

    #take a random sample from the posterior distribution
    pred_lik_pred_dens <- matrix(NA, nrow = nrow(pred), ncol = length(ind_samples))
    pred_lik_auc <- matrix(NA, ncol = length(ind_samples))
    
    for (m in 1:ncol(pred)) {
      pred_lik_pred_dens[,m] <- dpois(model_fits$training_data$response,
                                      pred[,m]*model_fits$training_data$weights,
                                      log = TRUE)
      
      occ_prob <- 1-exp(-exp(log(pred[,m]) + log(model_fits$training_data$weights)))
      suppressMessages(
        myRoc <- pROC::roc(response = model_fits$training_data$response, predictor = occ_prob))          
      pred_lik_auc[m] <- myRoc$auc
    }
    
    pred_dens <- apply(pred_lik_pred_dens, 1, function(x) -log(ncol(pred_lik_pred_dens)) + log(sum(exp(x))))
    pred_dens_int[i] <- sum(pred_dens)
    AUC_int[i] <- mean(pred_lik_auc)
  }
  
  thres_all <- list(TSS_thres, TP05_thres, AUC_int, pred_dens_int)
  names(thres_all) <- c('TSS', 'TP05', 'AUC_int', 'pred_dens')
  return(thres_all)
}