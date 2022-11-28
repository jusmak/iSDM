## Functions for computing predictive accuracy
#auc
return_auc <- function(response, pred_temp) {
  my_roc <- suppressMessages(pROC::roc(response = response, predictor = pred_temp))
  return(my_roc$auc)
}

#log-intensity value which maximizes sum of specificity and sensitivity
return_thres <- function(response, pred_temp) {
  my_roc <- suppressMessages(pROC::roc(response = response, predictor = pred_temp))
  return(my_roc$thresholds[which.max(my_roc$sensitivities+my_roc$specificities)[1]])
}

#kappa value
return_kappa <- function(response, pred_temp) {
  ck <- kappa2(cbind(response,pred_temp))
  return(ck$value)
}

#tss and other values
return_tss <- function(response, threshold, pred_temp) {
  my_roc <- suppressMessages(pROC::roc(response = response, predictor = pred_temp))
  ind_thres <- which.min(abs(my_roc$thresholds-threshold)) 
  tss <- max(my_roc$sensitivities[ind_thres]+my_roc$specificities[ind_thres]-1)
  spec <- my_roc$specificities[ind_thres]
  sens <- my_roc$sensitivities[ind_thres]
  return(list(tss, spec, sens))
}

#function for computing all predictive metrics
pred_metrics <- function(pred_pa, true_pa, pred_thres, true_thres) {
  pred_prob <- 1-exp(-exp(pred_pa))

  #auc
  auc_all <- mean(apply(pred_prob, 2, function(x) return_auc(true_pa, x)))
  
  #add a small jitter to allow computing likelihoods
  pred_prob_temp <- pred_prob
  pred_prob_temp[pred_prob_temp==0] <- pred_prob_temp[pred_prob_temp==0] + 1e-6
  pred_prob_temp[pred_prob_temp==1] <- pred_prob_temp[pred_prob_temp==1] - 1e-6
  
  log_lik <- apply(pred_prob_temp, 2, function(x) dbinom(true_pa, 1, x, log = T))
  log_lik_all <- sum(log(apply(log_lik, 1, function(x) mean(exp(x)))))
  tjur_r_all <- mean(apply(pred_prob, 2, function(x) mean(x[true_pa==1]) - mean(x[true_pa==0])))
  
  #threshold dependent metrics: kappa, tss, sensitivity, specificity
  thres_temp <- apply(pred_thres, 2, function(x) return_thres(response=true_thres, pred_temp=x))

  pred_cat <- matrix(NA, nrow = length(true_pa), ncol = ncol(pred_pa))
  for (l in 1:ncol(pred_thres)) {
    pred_cat[,l] <- ifelse(pred_pa[,l] < thres_temp[l], 0, 1)
  }  
  
  #kappa
  kappa_all_th <- mean(apply(pred_cat, 2, function(x) return_kappa(true_pa, x)))

  #tss, sens, spec
  tss_temp  <- 0
  spec_temp  <- 0
  sens_temp  <- 0
  for (k in 1:ncol(pred_prob)) {
    temp_2 <- return_tss(true_pa, thres_temp[k], pred_pa[,k])
    tss_temp[k] <- temp_2[[1]]
    spec_temp[k] <- temp_2[[2]]
    sens_temp[k] <- temp_2[[3]]
  }
  tss_all_th <- mean(tss_temp)
  spec_all_th <- mean(spec_temp)
  sens_all_th <- mean(sens_temp)
  
  
  #threshold dependent metrics: kappa, tss, sensitivity, specificity (use probability of .5)
  thres_temp <- log(-log(.5))
  
  pred_cat <- matrix(NA, nrow = length(true_pa), ncol = ncol(pred_pa))
  for (l in 1:ncol(pred_thres)) {
    pred_cat[,l] <- ifelse(pred_pa[,l] < thres_temp, 0, 1)
  }  
  
  #kappa
  kappa_all_p <- mean(apply(pred_cat, 2, function(x) return_kappa(true_pa, x)))
  
  #tss, sens, spec
  tss_temp  <- 0
  spec_temp  <- 0
  sens_temp  <- 0
  for (k in 1:ncol(pred_prob)) {
    temp_2 <- return_tss(true_pa, thres_temp, pred_pa[,k])
    tss_temp[k] <- temp_2[[1]]
    spec_temp[k] <- temp_2[[2]]
    sens_temp[k] <- temp_2[[3]]
  }
  tss_all_p <- mean(tss_temp)
  spec_all_p <- mean(spec_temp)
  sens_all_p <- mean(sens_temp)
  
  return(list(auc_samp = mean(auc_all), kappa_samp_th = mean(kappa_all_th), kappa_samp_p = mean(kappa_all_p), tss_samp_th = mean(tss_all_th), tss_samp_p = mean(tss_all_p),
              spec_samp_th = mean(spec_all_th), spec_samp_p = mean(spec_all_p), sens_samp_th = mean(sens_all_th), sens_samp_p = mean(sens_all_p), log_lik_samp = sum(log_lik_all), tjur_r_samp = mean(tjur_r_all)))
}