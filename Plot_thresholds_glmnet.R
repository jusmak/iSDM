Plot_thresholds <- function(model_fits, training_data, offset_full){
  
  TSS_thres <- matrix(NA, nrow = length(offset_full))
  TP05_thres <- matrix(NA, nrow = length(offset_full))

  #add offsets
  for (i in 1:length(offset_full)) {
    glmnet_temp <- model_fits$glmnet_model[[i]]
    offset_temp <- log(offset_full[[i]][,1]) + log(offset_full[[i]][,2])
    
    x_temp <- model.matrix(~ poly(training_data$covariates, degree = 2, raw = TRUE))
    pred <- predict.glmnet(object = glmnet_temp, newx = x_temp,
                    newoffset = offset_temp)
    
    pred_mean <- exp(pred[,ncol(pred)] + log(training_data$weights))
    occ_prob_mean <- 1-exp(-(pred_mean))
    roc_temp <- pROC::roc(response = training_data$response, predictor = occ_prob_mean)
    TSS_thres[i] <- pred[which.max(roc_temp$sensitivities+roc_temp$specificities-1),ncol(pred)]
    TP05_thres[i] <- pred[which(roc_temp$sensitivities<=.95)[1],ncol(pred)]
  }
  
  thres_all <- list(TSS_thres, TP05_thres)
  names(thres_all) <- c('TSS', 'TP05')
  return(thres_all)
}