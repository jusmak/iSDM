Plot_thresholds <- function(model_fits, training_data, offset_full){
  
  TSS_thres <- matrix(NA, nrow = length(offset_full))
  TP05_thres <- matrix(NA, nrow = length(offset_full))

  #make a prediction in training locations
  pred_matrix <- cbind(training_data$covariates, training_data$covariates^2)
  
  #add offsets
  for (i in 1:length(offset_full)) {
    coeff <- model_fits$fit_hb_comb_all[[i]]$coefficient[-2]
    pred <- matrix(rep(coeff[1], nrow(pred_matrix)), nrow = nrow(pred_matrix), byrow = TRUE) + pred_matrix %*% coeff[-1]
    
    offset_temp <- log(offset_full[[i]][,1]) + log(offset_full[[i]][,2])
    pred_mean <- exp(pred + offset_temp + log(training_data$weights))
    occ_prob_mean <- 1-exp(-(pred_mean))
    suppressMessages(
    roc_temp <- pROC::roc(response = training_data$response, predictor = array(occ_prob_mean)))
    TSS_thres[i] <- log(pred_mean)[which.max(roc_temp$sensitivities+roc_temp$specificities-1)]
    TP05_thres[i] <- log(pred_mean)[which(roc_temp$sensitivities<=.95)[1]]
  }
  
  thres_all <- list(TSS_thres, TP05_thres)
  names(thres_all) <- c('TSS', 'TP05')
  return(thres_all)
}

Plot_thresholds_test_1 <- function(model_fits, training_data, offset_full){
  
  TSS_thres <- matrix(NA, nrow = length(offset_full))
  TP05_thres <- matrix(NA, nrow = length(offset_full))
  
  #make a prediction in training locations
  pred_matrix <- cbind(training_data$covariates, training_data$covariates^2)
  
  for (i in 1:length(model_fits$fit_hb_comb_all)) {
    coeff <- model_fits$fit_hb_comb_all[[i]]$coefficient
    pred <- matrix(rep(coeff[1] + coeff[2], nrow(pred_matrix)), nrow = nrow(pred_matrix), byrow = TRUE) + pred_matrix %*% coeff[c(-1,-2)]
    
    pred_mean <- exp(pred + log(training_data$weights))
    occ_prob_mean <- 1-exp(-(pred_mean))
    suppressMessages(
      roc_temp <- pROC::roc(response = training_data$response, predictor = array(occ_prob_mean)))
    TSS_thres[i] <- log(pred_mean)[which.max(roc_temp$sensitivities+roc_temp$specificities-1)]
    TP05_thres[i] <- log(pred_mean)[which(roc_temp$sensitivities<=.95)[1]]
  }
  
  thres_all <- list(TSS_thres, TP05_thres)
  names(thres_all) <- c('TSS', 'TP05')
  return(thres_all)
}

Plot_thresholds_test_2 <- function(model_fits, training_data, offset_full){
  
  TSS_thres <- matrix(NA, nrow = length(offset_full), ncol = length(model_fits$fit_hb_comb_all))
  TP05_thres <- matrix(NA, nrow = length(offset_full), ncol = length(model_fits$fit_hb_comb_all))
  
  #make a prediction in training locations
  pred_matrix <- cbind(training_data$covariates, training_data$covariates^2)
  
  for (i in 1:length(model_fits$fit_hb_comb_all)) {
    for (j in 1:length(offset_full)) {
      coeff <- model_fits$fit_hb_comb_all[[i]][[j]]$coefficient
      pred <- matrix(rep(coeff[1] + coeff[2], nrow(pred_matrix)), nrow = nrow(pred_matrix), byrow = TRUE) + pred_matrix %*% coeff[c(-1,-2)]
      
      pred_mean <- exp(pred + log(training_data$weights) + log(offset_full[[j]][,1]) + log(offset_full[[j]][,2]))
      occ_prob_mean <- 1-exp(-(pred_mean))
      suppressMessages(
        roc_temp <- pROC::roc(response = training_data$response, predictor = array(occ_prob_mean)))
      TSS_thres[j,i] <- log(pred_mean)[which.max(roc_temp$sensitivities+roc_temp$specificities-1)]
      TP05_thres[j,i] <- log(pred_mean)[which(roc_temp$sensitivities<=.95)[1]]
    }
  }
  
  thres_all <- list(TSS_thres, TP05_thres)
  names(thres_all) <- c('TSS', 'TP05')
  return(thres_all)
}