Predictive_tests <- function(wd, model_fits, species, training_data, offset_full, pred_PA_data, proj_raster,
                             HB_inf, glmnet_inf, ppml_inf, thinning, target_n_obs, m_category,n_samp) {

  samp_temp <- extract(model_fits$HB_model[[1]])
  samp_temp <- sample(1:length(samp_temp$alpha),n_samp,replace=FALSE)
  thin_mark <- ifelse(thinning, "thin", "not_thin")  

  if(file.exists(paste(wd, '/Model_fits/', species, '_validation_', thin_mark, "_quad_n_",
                       target_n_obs, '_m_cat_', m_category, '.RData', sep = ''))) {
    
    load(paste(wd, '/Model_fits/', species, '_validation_', thin_mark, "_quad_n_",
               target_n_obs, '_m_cat_', m_category, '.RData', sep = ''))
  } else {

    #update here the real ind_fold
    #get occurrence points
    occ_points <- training_data$coordinates[training_data$response==1,]
    #transform the coordinates to the original scale
    occ_points <- cbind(occ_points[,1]*1000+training_data$min_coordinates[1],
                        occ_points[,2]*1000+training_data$min_coordinates[2])
    
    #create a data frame for species names
    sp_name <- data.frame(rep(species, nrow(occ_points)))
    colnames(sp_name) <- 'sp'
    
    #transform occurrence points into a spatialpoints data frame
    occ_sp_points <- SpatialPointsDataFrame(coords = occ_points,
                                            data = sp_name,
                                            proj4string = crs(proj_raster))
    
    #do a clustering-based division to folds
    source(paste(wd, "R_code/spatialCV.r", sep = '/'))
    sp_stats <- list(species)
    names(sp_stats) <- 'species'
    n_fold <- spatialStratify(occ_sp_points, sp_stats, nfolds = -1, nsubclusters = -1)
    ind_fold_temp <- n_fold@data$folds
    #create a list of folds
    ind_fold <- list()
    for (i in 1:length(unique(ind_fold_temp))) {
      ind_fold[[i]] <- which(ind_fold_temp == i)
    }
  
    ##HB
    if (HB_inf) {
      ##fit to data
      source(paste(wd, "R_code/Inference_HB_interp.R", sep = '/'))
      HB_train_valid <- Inference_HB_interp(wd, model_fits, training_data, offset_full, samp_temp)
      
      ##5-fold-cv
      source(paste(wd, "R_code/Inference_HB_cv_x_fold.R", sep = '/'))
      HB_cv_valid <- Inference_HB_cv_x_fold(wd, model_fits, training_data, offset_full, ind_fold)
    
      ##predictive density in the independent validation points
      source(paste(wd, "R_code/Inference_HB_ind_validation.R", sep = '/'))
      HB_ind_valid <- Inference_HB_ind_validation(wd, model_fits, pred_PA_data, samp_temp)
    } else {
      HB_train_valid <- list()
      HB_cv_valid <- list()
      HB_ind_valid <- list()
    }
  
    ##glmnet
    if (glmnet_inf) {
      ##fit to data
      source(paste(wd, "R_code/Inference_glmnet_interp.R", sep = '/'))
      glmnet_train_valid <- Inference_glmnet_interp(wd, model_fits, training_data, offset_full)
      
      ##5-fold-cv
      #fill here the folds in foldid as a vector indexing the fold where an observation belongs to
      source(paste(wd, "R_code/Inference_glmnet_cv_x_fold.R", sep = '/'))
      glmnet_cv_valid <- Inference_glmnet_cv_x_fold(wd, training_data, offset_full, ind_fold)
      
      ##predictive density in the independent validation points
      source(paste(wd, "R_code/Inference_glmnet_ind_validation.R", sep = '/'))
      glmnet_ind_valid <- Inference_glmnet_ind_validation(wd, model_fits, pred_PA_data)
    } else {
      glmnet_train_valid <- list()
      glmnet_cv_valid <- list()
      glmnet_ind_valid <- list()
    }
    
  
    ##ppmlasso
    if (ppml_inf) {
      source(paste(wd, "R_code/Inference_ppml_interp.R", sep = '/'))
      ppmlasso_train_valid <- Inference_ppml_interp(wd, model_fits, training_data, offset_full)
      
      ##5-fold-cv
      #fill here the folds in foldid as a vector indexing the fold where an observation belongs to
      source(paste(wd, "R_code/Inference_ppml_cv_x_fold.R", sep = '/'))
      ppmlasso_cv_valid <- Inference_ppml_cv_x_fold(wd, training_data, offset_full, ind_fold)
      
      ##predictive density in the independent validation points
      source(paste(wd, "R_code/Inference_ppml_ind_validation.R", sep = '/'))
      ppmlasso_ind_valid <- Inference_ppml_ind_validation(wd, model_fits, pred_PA_data)
    } else {
      ppmlasso_train_valid <- list()
      ppmlasso_cv_valid <- list()
      ppmlasso_ind_valid <- list()
    }
    
    predictive_validation <- list(HB_train_valid, HB_cv_valid, HB_ind_valid, glmnet_train_valid, glmnet_cv_valid, glmnet_ind_valid,
                                  ppmlasso_train_valid, ppmlasso_cv_valid, ppmlasso_ind_valid)
    names(predictive_validation) <- c('HB_train_valid', 'HB_cv_valid', 'HB_ind_valid', 'glmnet_train_valid', 'glmnet_cv_valid', 'glmnet_ind_valid',
                                      'ppmlasso_train_valid', 'ppmlasso_cv_valid', 'ppmlasso_ind_valid')
    
    save(predictive_validation, file = paste(wd, '/Model_fits/', species, '_validation_', thin_mark, "_quad_n_",
                                   target_n_obs, '_m_cat_', m_category, '.RData', sep = ''))  
  }
  return(predictive_validation)
  }


