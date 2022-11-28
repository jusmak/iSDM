RunInference = function(wd, species, data_inla, thinning, thinning_PA, weights_area, target_n_obs, model_v) {

  thin_mark = ifelse(thinning, 'thin', 'not_thin')
  thin_PA_mark =  ifelse(thinning_PA, '_thin_PA', '_not_thin_PA')
  weights_mark = ifelse(weights_area, 'wa', 'not_wa')
  
  if(file.exists(paste0('Model_fits_v', model_v, '/INLA_non_spat_', species, '_fits_', thin_mark, thin_PA_mark, '_quad_n_',
                       target_n_obs, '_', weights_mark, '_inla.RData'))) {
    
    load(paste0('Model_fits_v', model_v, '/INLA_non_spat_', species, '_fits_', thin_mark, thin_PA_mark, '_quad_n_',
               target_n_obs, '_', weights_mark, '_inla.RData'))
    
  } else {
    
    model_titles = c('non_spat_po', 'non_spat_pa', 'non_spat_comb')
    model_fits = vector(mode = 'list', length = length(model_titles))
    stk_coll = vector(mode = 'list', length = length(model_titles))
    names(model_fits) = model_titles
    names(stk_coll) = model_titles
    
    # set priors for the fixed effects
    

    #Each model building section contains following steps:
    # 1. set spde
    # 2. set data stack for PO data
    # 3. set data stack for PA data
    # 4. combine data stacks and define model structure
    # 5. fit model
    # 6. store data stack and model fit
    
    ##PO non-spatial model
    # 1.
    # ---
    
    # 2.
    stk_po = inla.stack(data = list(y=data_inla$po_response, e=data_inla$po_weights), A = list(1), tag = 'po',
                         effects = list(list(data.frame(Constant = rep(1, length(data_inla$po_response)), data_inla$po_covariates))))

    # 3.
    # ---
    
    # 4.
    f1 = as.formula(paste0('y ~ -1 + Constant + ', paste0(colnames(data_inla$po_covariates), collapse = ' + ')))
    
    # 5.
    ft.inla = 0
    ft.inla = try(inla(f1, family = 'poisson', data = inla.stack.data(stk_po),
                   control.fixed = list(prec = list(default=0.1), prec.intercept = list(default=0.1)),
                   control.family = list(list(link = 'log')),
                   control.predictor = list(A = inla.stack.A(stk_po), compute = TRUE),
                   E = inla.stack.data(stk_po)$e, control.compute = list(dic = TRUE),
                   verbose = TRUE))
    
    # 6.
    model_fits[[1]] = ft.inla
    stk_coll[[1]] = stk_po
    
    ##PA non-spatial model
    # 1.
    # ---
    
    # 2.
    # ---
    
    # 3.
    stk_pa = inla.stack(data = list(y=data_inla$pa_response, Ntrials=rep(1, length(data_inla$pa_response))), A = list(1), tag = 'pa',
                         effects = list(list(data.frame(Constant = rep(1, length(data_inla$pa_response)), data_inla$pa_covariates))))
    
    # 4.
    f1 = as.formula(paste0('y ~ -1 + Constant + ', paste0(colnames(data_inla$po_covariates), collapse = ' + ')))
    
    # 5.
    ft.inla = 0
    ft.inla = try(inla(f1, family = 'binomial', data = inla.stack.data(stk_pa),
                   control.fixed = list(prec = list(default=0.1), prec.intercept = list(default=0.1)),
                   control.family = list(list(link = 'cloglog')),
                   control.predictor = list(A = inla.stack.A(stk_pa), compute = TRUE),
                   Ntrials = inla.stack.data(stk_pa)$Ntrials, control.compute = list(dic = TRUE),
                   verbose = TRUE))
    
    # 6.
    model_fits[[2]] = ft.inla
    stk_coll[[2]] = stk_pa
    
    ##Combined non-spatial model
    # 1.
    # ---
    
    # 2.
    stk_po = inla.stack(data = list(y=cbind(data_inla$po_response, NA), e=data_inla$po_weights), A = list(1), tag = 'po',
                         effects = list(list(data.frame(ConstantB = rep(1, length(data_inla$po_response)), data_inla$po_covariates))))
    
    # 3.
    stk_pa = inla.stack(data = list(y=cbind(NA, data_inla$pa_response), Ntrials=rep(1, length(data_inla$pa_response))), A = list(1), tag = 'pa',
                         effects = list(list(data.frame(ConstantA = rep(1, length(data_inla$pa_response)), data_inla$pa_covariates))))
    
    # 4.
    stk_all = inla.stack(stk_po, stk_pa)
    f1 = as.formula(paste0('y ~ -1 + ConstantA + ConstantB + ', paste0(colnames(data_inla$po_covariates), collapse = ' + ')))
    
    # 5.
    ft.inla = 0
    ft.inla = try(inla(f1, family = c('poisson', 'binomial'), data = inla.stack.data(stk_all),
                   control.fixed = list(prec = list(default=0.1), prec.intercept = list(default=0.1)),
                   control.family = list(list(link = 'log'), list(link = 'cloglog')),
                   control.predictor = list(A = inla.stack.A(stk_all), link = 1, compute = TRUE),
                   Ntrials = inla.stack.data(stk_all)$Ntrials, E = inla.stack.data(stk_all)$e,
                   control.compute = list(dic = TRUE),
                   verbose = TRUE))
    
    # 6.
    model_fits[[3]] = ft.inla
    stk_coll[[3]] = stk_all
    
    output = list(model_fits = model_fits, stk_coll = stk_coll)
    
    save(output, file = paste0('Model_fits_v', model_v, '/INLA_non_spat_', species, '_fits_', thin_mark, thin_PA_mark, '_quad_n_',
                                  target_n_obs, '_', weights_mark, '_inla.RData'))
    
  } 
    
  return(output)
}