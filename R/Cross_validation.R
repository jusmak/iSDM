Validation_cv = function(wd, species, data, data_inla, thinning, thinning_PA, weights_area, target_n_obs, n_samples, model_v) {
  
  thin_mark = ifelse(thinning, 'thin', 'not_thin')
  thin_PA_mark =  ifelse(thinning_PA, '_thin_PA', '_not_thin_PA')
  weights_mark = ifelse(weights_area, 'wa', 'not_wa')
  
  if(file.exists(paste0('Model_fits_v', model_v, '/INLA_', species, '_validation_', thin_mark, thin_PA_mark, '_quad_n_',
                       target_n_obs, '_', weights_mark, '_inla.RData'))) {
    
    load(paste0('Model_fits_v', model_v, '/INLA_', species, '_validation_', thin_mark, thin_PA_mark, '_quad_n_',
               target_n_obs, '_', weights_mark, '_inla.RData'))
    
  } else {
    
    #read validation functions
    source('R/Validation_functions.R')
    
    #for INLA, create a mesh for approximating the Gaussian random field
    output_validation_cv = list()

    model_titles = c('non_spat_po', 'non_spat_pa', 'non_spat_comb', 'spat_po', 'spat_pa', 'spat_comb', 'spat_comb_pref',
                      'spat_po_rsr', 'spat_pa_rsr', 'spat_comb_rsr')
    
    #loop over prediction scenarios
    for (i in 1:length(data_inla$pa_testing_folds)) {
      
      model_valid = list()
      
      #loop over models
      for (j in 1:length(model_titles)) {

        if (j == 1) {
          load(file = paste0('Model_fits_v', model_v, '/INLA_non_spat_', species, '_fits_CV_', thin_mark, thin_PA_mark, '_quad_n_',
                                       target_n_obs, '_', weights_mark, '_inla.RData'))
        } else if (j == 4) {
          load(file = paste0('Model_fits_v', model_v, '/INLA_spat_', species, '_fits_CV_', thin_mark, thin_PA_mark, '_quad_n_',
                            target_n_obs, '_', weights_mark, '_inla.RData'))
        } else if (j == 7) {
          load(file = paste0('Model_fits_v', model_v, '/INLA_spat_samp_', species, '_fits_CV_', thin_mark, thin_PA_mark, '_quad_n_',
                            target_n_obs, '_', weights_mark, '_inla.RData'))
        } else if (j == 8) {
          load(file = paste0('Model_fits_v', model_v, '/INLA_spat_rsr_', species, '_fits_CV_', thin_mark, thin_PA_mark, '_quad_n_',
                            target_n_obs, '_', weights_mark, '_inla.RData'))
        }
        
        model_ind = c(1:3,1:3,1,1:3)[j]
        
        if (length(output_cv[[i]]$model_fits[[model_ind]]) > 1) {
          
          #find the prediction cells
          ind_pred = inla.stack.index(output_cv[[i]]$stk_coll[[model_ind]], 'pa_test')$data
          #find the training cells
          ind_thres = inla.stack.index(output_cv[[i]]$stk_coll[[model_ind]], 'pa')$data
          #combine indices
          ind_all = c(ind_pred, ind_thres)
          #collect their marginal density functions
          fitted_dist = output_cv[[i]]$model_fits[[model_ind]]$marginals.linear.predictor[ind_all]
          #take n_samples samples from each density function
          temp_list = lapply(fitted_dist, function(x) inla.rmarginal(n_samples, x))
          pred_all = matrix(unlist(temp_list), ncol=length(temp_list[[1]]), byrow=T)
          #split to prediction and thresholding sets
          pred_pa = pred_all[1:length(ind_pred),]
          pred_thres = pred_all[(length(ind_pred)+1):length(ind_all),]
          
          valid_cv = pred_metrics(pred_pa=pred_pa, true_pa=data_inla$pa_response[data_inla$pa_testing_folds[[i]]],
                                   pred_thres=pred_thres, true_thres=data_inla$pa_response[-data_inla$pa_testing_folds[[i]]])
          
          model_valid[[j]] = valid_cv
          
        } else {
          
          model_valid[[j]] = 0
          
        }
      }
      
      names(model_valid) = model_titles
      output_validation_cv[[i]] = model_valid
      
    }
    names(output_validation_cv) = names(data_inla$pa_testing_folds)
    save(output_validation_cv, file = paste0('Model_fits_v', model_v, '/INLA_', species, '_validation_', thin_mark, thin_PA_mark, '_quad_n_',
                                            target_n_obs, '_', weights_mark, '_inla.RData'))
  }
}



