RunInference_v1 <- function(wd, species, data, thinning, thinning_PA, weights_area, target_n_obs, m_category) {

  thin_mark <- ifelse(thinning, "thin", "not_thin")
  thin_PA_mark <-  ifelse(thinning_PA, "_thin_PA", "_not_thin_PA")
  weights_mark <- ifelse(weights_area, "wa", "not_wa")
  if(file.exists(paste(wd, "/Model_fits_test_1/", species, '_fits_', thin_mark, thin_PA_mark, "_quad_n_",
                       target_n_obs, "_", weights_mark, '_m_cat_', m_category, '.RData', sep = ''))) {
    
    load(paste(wd, "/Model_fits_test_1/", species, '_fits_', thin_mark, thin_PA_mark, "_quad_n_",
               target_n_obs, "_", weights_mark, '_m_cat_', m_category, '.RData', sep = ''))
  } else {
    
    fit_hb_comb_all <- list()
    source(paste(wd, "R_code/Inference_HB_comb_test_1.R", sep = '/'))
    
    n_pa <- length(data$inventory_data$response)
    n_po <- sum(data$training_data$weights)
    
    #change the relative importance of likelihood of PO-data compared to PA-data
    weights_temp_PO <- c(0, .33, .5, .66, 1)
    weights_temp_PA <- rev(weights_temp_PO)
    
    #set the scaling scheme
    if (m_category==1) {
      scale_temp <- 1
    }
    if (m_category==2) {
      scale_temp <- n_pa/n_po
    }
    if (m_category==3) {
      scale_temp <- log(n_pa)/log(n_po)
    }
    
    
    #conduct Empirical Bayesian inference
    for (j in 1:length(weights_temp_PA)) {
      data_temp <- data
      
      #scale weights of the PO-data
      data_temp$training_data <- append(data_temp$training_data,
                                   list(rep(scale_temp*weights_temp_PO[j],length(data_temp$training_data$response))))
      names(data_temp$training_data)[length(data_temp$training_data)] <- "weights_lik"
      
      weights_list <- rep(weights_temp_PA[j], length(data_temp$inventory_data$response))
      data_temp$inventory_data <- append(data_temp$inventory_data,
                                         list(rep(weights_temp_PA[j], n_pa)))
      names(data_temp$inventory_data)[length(data_temp$inventory_data)] <- "weights"
      
      fit_hb_weights <- list()
      
      for (k in 1:length(data_temp$offset_full)) {
        fit_hb <- Inference_HB_comb(wd, data_temp$training_data, data_temp$inventory_data,
                                    data_temp$offset_full[[k]], data_temp$inventory_data$offset[[k]])
        fit_hb_weights[[k]] <- fit_hb
      }
        
      fit_hb_comb_all[[j]] <- fit_hb_weights
    }
    
    model_fits <- list(fit_hb_comb_all)
    names(model_fits) <- c('fit_hb_comb_all')
    
    save(model_fits, file = paste(wd, "/Model_fits_test_1/", species, '_fits_', thin_mark, thin_PA_mark, "_quad_n_",
                                  target_n_obs, "_", weights_mark, '_m_cat_', m_category, '.RData', sep = ''))
  }
  return(model_fits)
}