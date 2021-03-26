RunInference_v1 <- function(wd, species, data, thinning, thinning_PA, weights_area, target_n_obs, m_category) {

  thin_mark <- ifelse(thinning, "thin", "not_thin")
  thin_PA_mark <-  ifelse(thinning_PA, "_thin_PA", "_not_thin_PA")
  weights_mark <- ifelse(weights_area, "wa", "not_wa")
  if(file.exists(paste(wd, "/Model_fits_test_2/", species, '_fits_', thin_mark, thin_PA_mark, "_quad_n_",
                       target_n_obs, "_", weights_mark, '_m_cat_', m_category, '.RData', sep = ''))) {
    
    load(paste(wd, "/Model_fits_test_2/", species, '_fits_', thin_mark, thin_PA_mark, "_quad_n_",
               target_n_obs, "_", weights_mark, '_m_cat_', m_category, '.RData', sep = ''))
  } else {
    
    fit_hb_comb_all <- list()
    fit_hb_single_all <- list()
    
    source(paste(wd, "R_code/Inference_HB_comb_test_2.R", sep = '/'))
    
    #create indexes for the background points
    n_temp <- c(1, 50, 100, 200, 500, 1000, 5000, 10000, sum(data$training_data$response==0))

    #conduct Empirical Bayesian inference
    ind_pres <- which(data$training_data$response==1)
    
    for (j in 1:length(n_temp)) {
      data_temp <- data
      
      set.seed(j)
      ind_temp <- c(ind_pres, sort(sample(which(data_temp$training_data$response==0), n_temp[j], replace = FALSE)))
      data_temp$training_data$coordinates <- data_temp$training_data$coordinates[ind_temp,]
      data_temp$training_data$covariates <- data_temp$training_data$covariates[ind_temp,]
      data_temp$training_data$response <- data_temp$training_data$response[ind_temp]
      data_temp$training_data$weights <- data_temp$training_data$weights[ind_temp]
      data_temp$offset_full <- lapply(data_temp$offset_full, function(x) x[ind_temp,])
      
      ind_background <- which(data_temp$training_data$response==0 & data_temp$training_data$weights==.5)
      ind_pres_sub <- which(data_temp$training_data$response==1)
#      if (length(ind_background)==0) {
#        data_temp$training_data$weights <- rep(1, length(data_temp$training_data$weights))
#      } else {
#        ind_2 <- apply(as.matrix(ind_background), 1, function(x) which.min((data_temp$training_data$coordinates[x,1] - 
#                                                                  data_temp$training_data$coordinates[ind_pres_sub,1])^2 +
#                                                                (data_temp$training_data$coordinates[x,2] - 
#                                                                   data_temp$training_data$coordinates[ind_pres_sub,2])^2))
#        data_temp$training_data$weights[ind_pres_sub[-ind_2]] <- rep(1, length(ind_pres_sub[-ind_2]))
#      }
      data_temp$training_data$weights <- rep(sum(data$training_data$weights[data$training_data$response==0])/n_temp[j], 
                                             length(data_temp$training_data$weights))
      data_temp$training_data$weights[ind_pres_sub] = rep(1e-6,length(ind_pres_sub))

      fit_hb <- Inference_HB_comb(wd, data_temp$training_data, data_temp$inventory_data,
                                  data_temp$offset_full[[1]], data_temp$inventory_data$offset[[1]])
      fit_hb_comb_all[[j]] <- fit_hb
      
      fit_hb <- Inference_HB_single(wd, data_temp$training_data, data_temp$inventory_data,
                                  data_temp$offset_full[[1]], data_temp$inventory_data$offset[[1]])
      fit_hb_single_all[[j]] <- fit_hb
    }
    
    model_fits <- list(fit_hb_comb_all, fit_hb_single_all)
    names(model_fits) <- c('fit_hb_comb_all', 'fit_hb_single_all')
    
    save(model_fits, file = paste(wd, "/Model_fits_test_2/", species, '_fits_', thin_mark, thin_PA_mark, "_quad_n_",
                                  target_n_obs, "_", weights_mark, '_m_cat_', m_category, '.RData', sep = ''))
  }
  return(model_fits)
}