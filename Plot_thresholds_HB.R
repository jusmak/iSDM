Plot_thresholds <- function(model_fits, training_data, offset_full, ind_samples){
  
  TSS_thres <- matrix(NA, nrow = length(offset_full))
  TP05_thres <- matrix(NA, nrow = length(offset_full))

  #make a prediction in training locations
  pred_matrix <- cbind(training_data$covariates, training_data$covariates^2)
  
  #add offsets
  for (i in 1:length(offset_full)) {
    beta_samp <- extract(model_fits$HB_model[[i]])
    pred <- matrix(rep(beta_samp$alpha[ind_samples], nrow(pred_matrix)), nrow = nrow(pred_matrix), byrow = TRUE) + pred_matrix %*% t(beta_samp$beta[ind_samples,])
    
    offset_temp <- log(offset_full[[i]][,1]) + log(offset_full[[i]][,2])
    pred <- exp(pred + offset_temp + log(training_data$weights))
    pred_mean <- apply(pred, 1, mean)
    occ_prob_mean <- 1-exp(-(pred_mean))
    roc_temp <- pROC::roc(response = training_data$response, predictor = occ_prob_mean)
    TSS_thres[i] <- log(pred)[which.max(roc_temp$sensitivities+roc_temp$specificities-1)]
    TP05_thres[i] <- log(pred)[which(roc_temp$sensitivities<=.95)[1]]
  }
  
  thres_all <- list(TSS_thres, TP05_thres)
  names(thres_all) <- c('TSS', 'TP05')
  return(thres_all)
}