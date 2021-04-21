Predictive_tests <- function(wd, model_fits, species, training_data, offset_full, pred_PA_data, 
                             proj_raster, thinning, thinning_PA, weights_area, target_n_obs, m_category, n_samp) {

  thin_mark <- ifelse(thinning, "thin", "not_thin")
  thin_PA_mark <-  ifelse(thinning_PA, "_thin_PA", "_not_thin_PA")
  weights_mark <- ifelse(weights_area, "wa", "not_wa")
  if(file.exists(paste(wd, '/Model_fits_test_3/', species, '_validation_', thin_mark, thin_PA_mark, "_quad_n_",
                       target_n_obs, "_", weights_mark, '_m_cat_', m_category, '.RData', sep = ''))) {
    
    load(paste(wd, '/Model_fits_test_3/', species, '_validation_', thin_mark, thin_PA_mark, "_quad_n_",
               target_n_obs, "_", weights_mark, '_m_cat_', m_category, '.RData', sep = ''))
  } else {

    # output matrices
    auc <- matrix(NA, ncol = length(model_fits$fit_hb_comb_all), nrow = 1)
    tss <- matrix(NA, ncol = length(model_fits$fit_hb_comb_all), nrow = 1)
    spec <- matrix(NA, ncol = length(model_fits$fit_hb_comb_all), nrow = 1)
    sens <- matrix(NA, ncol = length(model_fits$fit_hb_comb_all), nrow = 1)
    kappa <- matrix(NA, ncol = length(model_fits$fit_hb_comb_all), nrow = 1)
    pred_dens <- matrix(NA, ncol = length(model_fits$fit_hb_comb_all), nrow = 1)
    tjur_r <- matrix(NA, ncol = length(model_fits$fit_hb_comb_all), nrow = 1)
    
    pred_matrix <- cbind(pred_PA_data$covariates, pred_PA_data$covariates^2)
    
    #loop over different models (offsets used for model fitting)
    for (i in 1:length(model_fits$fit_hb_comb_all)) {
      beta_samp <- model_fits$fit_hb_comb_all[[i]]$coefficient
      pred_validation <- matrix(rep(beta_samp[1] + beta_samp[2], nrow(pred_matrix)), nrow = nrow(pred_matrix), byrow = TRUE) + pred_matrix %*% beta_samp[3:length(beta_samp)]
      pred_prob <- 1-exp(-exp(pred_validation))
      
      suppressMessages(
        myRoc <- pROC::roc(response = pred_PA_data$response, predictor = as.vector(1-exp(-exp(pred_validation)))))
      auc[i] <- myRoc$auc
      pred_cat <- ifelse(pred_prob < .5, 0, 1)
      ck <- kappa2(cbind(pred_PA_data$response,pred_cat))
      kappa[i] <- ck$value
      ind_thres <- which.min(abs(myRoc$thresholds-.5)) 
      tss[i] <- max(myRoc$sensitivities[ind_thres]+myRoc$specificities[ind_thres]-1)
      spec[i] <- myRoc$specificities[ind_thres]
      sens[i] <- myRoc$sensitivities[ind_thres]
      pred_dens[i] <- sum(dbinom(pred_PA_data$response, 1,pred_prob, log = TRUE))
      tjur_r[i] <- mean(pred_prob[pred_PA_data$response==1]) - mean(pred_prob[pred_PA_data$response==0])
    }
    
    model_valid <- list(auc, tss, spec, sens, kappa, pred_dens, tjur_r)
    names(model_valid) <- c('auc', 'tss', 'specificity', 'sensitivity', 'kappa', 'pred_dens', 'tju_r')
    save(model_valid, file = paste(wd, '/Model_fits_test_3/', species, '_validation_', thin_mark, thin_PA_mark, "_quad_n_",
                                             target_n_obs, "_", weights_mark, '_m_cat_', m_category, '.RData', sep = ''))
  }
  return(model_valid)
  }


