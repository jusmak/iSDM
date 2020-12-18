RunInference_v1 <- function(wd, species, data, HB_inf, glmnet_inf, ppml_inf, lik, thinning, target_n_obs, m_category) {

  thin_mark <- ifelse(thinning, "thin", "not_thin")  
  if(file.exists(paste(wd, "/Model_fits/", species, '_fits_', thin_mark, "_quad_n_",
                       target_n_obs, '_m_cat_', m_category, '.RData', sep = ''))) {
    
    load(paste(wd, "/Model_fits/", species, '_fits_', thin_mark, "_quad_n_",
               target_n_obs, '_m_cat_', m_category, '.RData', sep = ''))
  } else {
      
    #conduct Hierarchical Bayesian inference for each different offset scenario
    fit_hb_all <- list()
    if (HB_inf) {
      source(paste(wd, "R_code/Inference_HB.R", sep = '/'))
      for (i in 1:4) {
        fit_hb <- Inference_HB(wd, data$training_data, data$offset_full[[i]])
        fit_hb_all[[i]] <- fit_hb
      }
    }
    
    #conduct frequentist regularized regression for each different offset scenario
    fit_glmnet_all <- list()
    if (glmnet_inf) {
      source(paste(wd, "R_code/Inference_glmnet.R", sep = '/'))
      for (i in 1:4) {
        fit_glmnet <- Inference_glmnet(wd, data$training_data, data$offset_full[[i]])
        fit_glmnet_all[[i]] <- fit_glmnet
      }
    }
    
    #conduct frequentist regularized regression with ppmlasso for each different offset scenario
    fit_ppml_all <- list()
    fit_ppml_comb_all <- list()
    if (ppml_inf) {
      source(paste(wd, "R_code/Inference_ppml.R", sep = '/'))
      for (i in 1) {
        fit_ppml <- Inference_ppml(wd, data$training_data, data$offset_full[[i]])
        fit_ppml_all[[i]] <- fit_ppml
      }
      if (lik_comb) {
        source(paste(wd, "R_code/Inference_ppml_comb.R", sep = '/'))
        for (i in 1) {
          fit_ppml <- Inference_ppml_comb(wd, data$training_data, data$inventory_data, data$offset_full[[i]])
          fit_ppml_comb_all[[i]] <- fit_ppml
        }
      }
    }
    
    model_fits <- list(fit_hb_all,fit_glmnet_all,fit_ppml_all, fit_ppml_comb_all)
    names(model_fits) <- c('HB_model', 'glmnet_model', 'ppml_model', 'ppml_comb_model')
    
    save(model_fits, file = paste(wd, "/Model_fits/", species, '_fits_', thin_mark, "_quad_n_",
                                  target_n_obs, '_m_cat_', m_category, '.RData', sep = ''))
  }
  return(model_fits)
}