Predictive_tests <- function(wd, model_fits, species, training_data, offset_full, pred_PA_data, 
                             proj_raster, thinning, thinning_PA, weights_area, target_n_obs, m_category, n_samp) {

  thin_mark <- ifelse(thinning, "thin", "not_thin")
  thin_PA_mark <-  ifelse(thinning_PA, "_thin_PA", "_not_thin_PA")
  weights_mark <- ifelse(weights_area, "wa", "not_wa")
  if(file.exists(paste(wd, '/Model_fits_test_1/', species, '_validation_', thin_mark, thin_PA_mark, "_quad_n_",
                       target_n_obs, "_", weights_mark, '_m_cat_', m_category, '.RData', sep = ''))) {
    
    load(paste(wd, '/Model_fits_test_1/', species, '_validation_', thin_mark, thin_PA_mark, "_quad_n_",
               target_n_obs, "_", weights_mark, '_m_cat_', m_category, '.RData', sep = ''))
  } else {

    
    # output matrices
    temp_col <- length(model_fits$fit_hb_comb_all) # different weights
    temp_row <- length(model_fits$fit_hb_comb_all[[1]]) # different offsets
    auc <- matrix(NA, ncol = temp_col, nrow = temp_row)
    tss <- matrix(NA, ncol = temp_col, nrow = temp_row)
    spec <- matrix(NA, ncol = temp_col, nrow = temp_row)
    sens <- matrix(NA, ncol = temp_col, nrow = temp_row)
    kappa <- matrix(NA, ncol = temp_col, nrow = temp_row)
    pred_dens <- matrix(NA, ncol = temp_col, nrow = temp_row)
    tjur_r <- matrix(NA, ncol = temp_col, nrow = temp_row)
    
    pred_matrix <- cbind(pred_PA_data$covariates, pred_PA_data$covariates^2)
    
    #loop over different models (offsets used for model fitting)
    for (j in 1:temp_row) {
      offset_temp <- log(pred_PA_data$offset[[j]][,1]) + log(pred_PA_data$offset[[j]][,2])
      
      #loop over different weights (used for model weighing datasets)
      for (i in 1:temp_col) {
        beta_samp <- model_fits$fit_hb_comb_all[[i]][[j]]$coefficient
        pred_validation <- matrix(rep(beta_samp[1] + beta_samp[2], nrow(pred_matrix)), nrow = nrow(pred_matrix), byrow = TRUE) + pred_matrix %*% beta_samp[3:length(beta_samp)]
        pred_prob <- 1-exp(-exp(pred_validation + offset_temp))
        
        suppressMessages(
          myRoc <- pROC::roc(response = pred_PA_data$response, predictor = as.vector(1-exp(-exp(pred_validation + offset_temp)))))
        auc[j,i] <- myRoc$auc
        pred_cat <- ifelse(pred_prob < .5, 0, 1)
        ck <- kappa2(cbind(pred_PA_data$response,pred_cat))
        kappa[j,i] <- ck$value
        ind_thres <- which.min(abs(myRoc$thresholds-.5)) 
        tss[j,i] <- max(myRoc$sensitivities[ind_thres]+myRoc$specificities[ind_thres]-1)
        spec[j,i] <- myRoc$specificities[ind_thres]
        sens[j,i] <- myRoc$sensitivities[ind_thres]
        pred_dens[j,i] <- sum(dbinom(pred_PA_data$response, 1,pred_prob, log = TRUE))
        tjur_r[j,i] <- mean(pred_prob[pred_PA_data$response==1]) - mean(pred_prob[pred_PA_data$response==0])
      }      
    }

    
    model_valid <- list(auc, tss, spec, sens, kappa, pred_dens, tjur_r)
    names(model_valid) <- c('auc', 'tss', 'specificity', 'sensitivity', 'kappa', 'pred_dens', 'tjur_r')
    save(model_valid, file = paste(wd, '/Model_fits_test_1/', species, '_validation_', thin_mark, thin_PA_mark, "_quad_n_",
                                             target_n_obs, "_", weights_mark, '_m_cat_', m_category, '.RData', sep = ''))
  }
  return(model_valid)
  }


