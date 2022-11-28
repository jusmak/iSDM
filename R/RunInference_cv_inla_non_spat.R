RunInference_cv = function(wd, species, data_inla, thinning, thinning_PA, weights_area, target_n_obs, model_v) {

  thin_mark = ifelse(thinning, 'thin', 'not_thin')
  thin_PA_mark =  ifelse(thinning_PA, '_thin_PA', '_not_thin_PA')
  weights_mark = ifelse(weights_area, 'wa', 'not_wa')
  
  if(file.exists(paste0('Model_fits_v', model_v, '/INLA_non_spat_', species, '_fits_CV_', thin_mark, thin_PA_mark, '_quad_n_',
                       target_n_obs, '_', weights_mark, '_inla.RData'))) {
    
    load(paste0('Model_fits_v', model_v, '/INLA_non_spat_', species, '_fits_CV_', thin_mark, thin_PA_mark, '_quad_n_',
               target_n_obs, '_', weights_mark, '_inla.RData'))
    
    #look for the validation set where there is no more model fit
    temp_i = matrix(NA, length(data_inla$pa_testing_folds), 1)
    for (k in 1:length(data_inla$pa_testing_folds)) {
      temp_i[k] = length(output_cv[[k]])
    }
    
    temp_i = which(temp_i==0)[1]
    
  } else {
    
    #store all fits and data stacks in a list
    output_cv = vector(mode = 'list', length = length(data_inla$pa_testing_folds))
    names(output_cv) = names(data_inla$pa_testing_folds)
    temp_i = 1
    
  }
  
  if (!is.na(temp_i)) {
    
    for (i in temp_i:length(data_inla$pa_testing_folds)) {
  
      model_titles = c('non_spat_po', 'non_spat_pa', 'non_spat_comb')
      model_fits = vector(mode = 'list', length = length(model_titles))
      stk_coll = vector(mode = 'list', length = length(model_titles))
      names(model_fits) = model_titles
      names(stk_coll) = model_titles
      
      #create a temporary list of data with pa data reduced to the training data set
      data_inla_temp = data_inla
      data_inla_temp$pa_coordinates = data_inla_temp$pa_coordinates[-data_inla_temp$pa_testing_folds[[i]],]
      data_inla_temp$pa_covariates = data_inla_temp$pa_covariates[-data_inla_temp$pa_testing_folds[[i]],]
      data_inla_temp$pa_response = data_inla_temp$pa_response[-data_inla_temp$pa_testing_folds[[i]]]
      
      #create the same lists but for testing data
      data_inla_temp$pa_coordinates_test = data_inla$pa_coordinates[data_inla_temp$pa_testing_folds[[i]],]
      data_inla_temp$pa_covariates_test = data_inla$pa_covariates[data_inla_temp$pa_testing_folds[[i]],]
      data_inla_temp$pa_response_test = data_inla$pa_response[data_inla_temp$pa_testing_folds[[i]]]
      
      #create training sets for the presence-only data
      data_inla_temp$po_coordinates = data_inla$po_coordinates[data_inla_temp$po_training_folds[[i]],]
      data_inla_temp$po_covariates = data_inla$po_covariates[data_inla_temp$po_training_folds[[i]],]
      data_inla_temp$po_response = data_inla$po_response[data_inla_temp$po_training_folds[[i]]]
      data_inla_temp$po_weights = data_inla$po_weights[data_inla_temp$po_training_folds[[i]]]
      
      if (length(data_inla_temp$po_response) == 0) {
        
      } else {
        
        #Each model building section contains following steps:
        # 1. set spde
        # 2. set data stack for PO data
        # 3. set data stack for PA data (train and test)
        # 4. combine data stacks and define model structure
        # 5. fit model
        # 6. store data stack and model fit
        
        ##PO non-spatial model
        # 1.
        # ---
        
        # 2.
        stk_po = inla.stack(data = list(y=cbind(data_inla_temp$po_response, NA), e=data_inla_temp$po_weights), A = list(1), tag = 'po',
                             effects = list(list(data.frame(Constant = rep(1, length(data_inla_temp$po_response)), data_inla_temp$po_covariates))))
        
        # 3.
        stk_pa_ref = inla.stack(data = list(y=matrix(NA, nrow = length(data_inla_temp$pa_response), ncol = 2), Ntrials=rep(1, length(data_inla_temp$pa_response))),
                                  A = list(1), tag = 'pa', 
                                  effects = list(list(data.frame(Constant = rep(1, length(data_inla_temp$pa_response)), data_inla_temp$pa_covariates))))
        
        stk_pa_test = inla.stack(data = list(y=matrix(NA, nrow = length(data_inla_temp$pa_response_test), ncol = 2), Ntrials=rep(1, length(data_inla_temp$pa_response_test))),
                                  A = list(1), tag = 'pa_test', 
                                  effects = list(list(data.frame(Constant = rep(1, length(data_inla_temp$pa_response_test)), data_inla_temp$pa_covariates_test))))
        
        # 4.
        stk_all = inla.stack(stk_po, stk_pa_test, stk_pa_ref)
        f1 = as.formula(paste0('y ~ -1 + Constant + ', paste0(colnames(data_inla_temp$po_covariates), collapse = ' + ')))
        
        # 5.
        ft.inla = 0
        ft.inla = try(inla(f1, family = c('poisson', 'binomial'), data = inla.stack.data(stk_all),
                       control.family = list(list(link = 'log'), list(link = 'cloglog')),
                       control.predictor = list(A = inla.stack.A(stk_all), link = 2, compute = TRUE),
                       E = inla.stack.data(stk_all)$e, Ntrials = inla.stack.data(stk_all)$Ntrials,
                       control.compute = list(dic = TRUE)))
        
        # 6.
        model_fits[[1]] = ft.inla
        stk_coll[[1]] = stk_all
        
        ##PA non-spatial model
        # 1.
        # ---
        
        # 2.
        # ---
        
        # 3.
        stk_pa = inla.stack(data = list(y=data_inla_temp$pa_response, Ntrials=rep(1, length(data_inla_temp$pa_response))), A = list(1), tag = 'pa',
                             effects = list(list(data.frame(Constant = rep(1, length(data_inla_temp$pa_response)), data_inla_temp$pa_covariates))))
        
        stk_pa_test = inla.stack(data = list(y=matrix(NA, nrow = length(data_inla_temp$pa_response_test), ncol = 1), Ntrials=rep(1, length(data_inla_temp$pa_response_test))),
                                  A = list(1), tag = 'pa_test', 
                                  effects = list(list(data.frame(Constant = rep(1, length(data_inla_temp$pa_response_test)), data_inla_temp$pa_covariates_test))))
        
        # 4.
        stk_all = inla.stack(stk_pa, stk_pa_test)
        f1 = as.formula(paste0('y ~ -1 + Constant + ', paste0(colnames(data_inla_temp$po_covariates), collapse = ' + ')))
        
        # 5.
        ft.inla = 0
        ft.inla = try(inla(f1, family = 'binomial', data = inla.stack.data(stk_all),
                       control.fixed = list(prec = list(default=0.1), prec.intercept = list(default=0.1)),
                       control.family = list(list(link = 'cloglog')),
                       control.predictor = list(A = inla.stack.A(stk_all), link = 1, compute = TRUE),
                       Ntrials = inla.stack.data(stk_all)$Ntrials, control.compute = list(dic = TRUE),
                       verbose = TRUE))
        
        # 6.
        model_fits[[2]] = ft.inla
        stk_coll[[2]] = stk_all
        
        ##Combined non-spatial model
        # 1.
        # ---
        
        # 2.
        stk_po = inla.stack(data = list(y=cbind(data_inla_temp$po_response, NA), e=data_inla_temp$po_weights), A = list(1), tag = 'po',
                             effects = list(list(data.frame(ConstantB = rep(1, length(data_inla_temp$po_response)), data_inla_temp$po_covariates))))
        
        # 3.
        stk_pa = inla.stack(data = list(y=cbind(NA, data_inla_temp$pa_response), Ntrials=rep(1, length(data_inla_temp$pa_response))), A = list(1), tag = 'pa',
                             effects = list(list(data.frame(ConstantA = rep(1, length(data_inla_temp$pa_response)), data_inla_temp$pa_covariates))))
  
        stk_pa_test = inla.stack(data = list(y=matrix(NA, nrow = length(data_inla_temp$pa_response_test), ncol = 2), Ntrials=rep(1, length(data_inla_temp$pa_response_test))),
                                  A = list(1), tag = 'pa_test', 
                                  effects = list(list(data.frame(ConstantA = rep(1, length(data_inla_temp$pa_response_test)), data_inla_temp$pa_covariates_test))))
        
        # 4.
        stk_all = inla.stack(stk_po, stk_pa, stk_pa_test)
        f1 = as.formula(paste0('y ~ -1 + ConstantA + ConstantB + ', paste0(colnames(data_inla_temp$po_covariates), collapse = ' + ')))
        
        # 5.
        ft.inla = 0
        ft.inla = try(inla(f1, family = c('poisson', 'binomial'), data = inla.stack.data(stk_all),
                       control.fixed = list(prec = list(default=0.1), prec.intercept = list(default=0.1)),
                       control.family = list(list(link = 'log'), list(link = 'cloglog')),
                       control.predictor = list(A = inla.stack.A(stk_all), link = 2, compute = TRUE),
                       Ntrials = inla.stack.data(stk_all)$Ntrials, E = inla.stack.data(stk_all)$e,
                       control.compute = list(dic = TRUE),
                       verbose = TRUE))
        
        # 6.
        model_fits[[3]] = ft.inla
        stk_coll[[3]] = stk_all
        
        output_cv[[i]] = list('model_fits' = model_fits, 'stk_coll' = stk_coll)
      }
      save(output_cv, file = paste0('Model_fits_v', model_v, '/INLA_non_spat_', species, '_fits_CV_', thin_mark, thin_PA_mark, '_quad_n_',
                                   target_n_obs, '_', weights_mark, '_inla.RData'))
    }
  }
  
  return(output_cv)
}
