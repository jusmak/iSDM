Plot_thresholds <- function(model_fits){
  
  TSS_thres <- matrix(NA, nrow = length(model_fits$offset))
  TP05_thres <- matrix(NA, nrow = length(model_fits$offset))
  AUC_int <- matrix(NA, nrow = length(model_fits$offset))
  pred_dens_int <- matrix(NA, nrow = length(model_fits$offset))
  
  #make a prediction in training locations
  pred_matrix <- cbind(model_fits$training_data$covariates, model_fits$training_data$covariates^2)
  
  #add offsets
  for (i in 1:length(model_fits$offset)) {
    glmnet_temp <- model_fits$glmnet_model[[i]]
    offset_temp <- log(model_fits$offset[[i]][,1]) - mean(log(model_fits$offset[[i]][,1])) +
      log(model_fits$offset[[i]][,2]) - mean(log(model_fits$offset[[i]][,2]))
    pred <- predict(object = glmnet_temp, newx = model_fits$training_data$covariates,
                    newoffset = offset_temp)
    
    pred_mean <- exp(pred[,ncol(pred)] + log(model_fits$training_data$weights))
    occ_prob_mean <- 1-exp(-(pred_mean))
    roc_temp <- pROC::roc(response = model_fits$training_data$response, predictor = occ_prob_mean)
    TSS_thres[i] <- pred[which.max(roc_temp$sensitivities+roc_temp$specificities-1),ncol(pred)]
    TP05_thres[i] <- pred[which(roc_temp$sensitivities<=.95)[1],ncol(pred)]

    #take a random sample from the posterior distribution
    pred_dens <- dpois(model_fits$training_data$response,
                                    exp(pred[,ncol(pred)]+log(model_fits$training_data$weights)),
                                    log = TRUE)
    suppressMessages(
      myRoc <- pROC::roc(response = model_fits$training_data$response, predictor = occ_prob_mean))          
    pred_dens_int[i] <- sum(pred_dens)
    AUC_int[i] <- myRoc$auc
  }
  
  thres_all <- list(TSS_thres, TP05_thres, AUC_int, pred_dens_int)
  names(thres_all) <- c('TSS', 'TP05', 'AUC_int', 'pred_dens')
  return(thres_all)
}