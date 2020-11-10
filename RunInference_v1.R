RunInference_v1 <- function(wd, species, observations, geographic_extent,thinning) {
  
  # get environmental data set
  source(paste(wd, "R_code/Get_env_data.R", sep = '/'))
  env_data <- Get_env_data(wd, geographic_extent)
  
  # get species observations
  source(paste(wd, "R_code/Get_species_data.R", sep = '/'))
  species_data <- Get_species_data(wd,species,env_data$min_coordinates,thinning)
  
  #sort data for inference
  source(paste(wd, "R_code/PP_training_data.R", sep = '/'))
  training_data <- PP_training_data(wd, env_data, species_data)
  
  #define offsets (expert range map + elevation)
  source(paste(wd, "R_code/Get_offsets.R", sep = '/'))
  offset_full <- Get_offsets(wd, training_data, env_data, species)

  #conduct Hierarchical Bayesian inference for each different offset scenario
  source(paste(wd, "R_code/Inference_HB.R", sep = '/'))
  fit_hb_all <- list()
  for (i in 1:4) {
    fit_hb <- Inference_HB(wd, training_data, offset_full[[i]])
    fit_hb_all[[i]] <- fit_hb
  }
  
  #conduct frequentist regularized regression
  source(paste(wd, "R_code/Inference_glmnet.R", sep = '/'))
  fit_glmnet_all <- list()
  for (i in 1:4) {
    fit_glmnet <- Inference_glmnet(wd, training_data, offset_full[[i]])
    fit_glmnet_all[[i]] <- fit_glmnet
  }
  
  model_fits <- list(fit_hb_all,fit_glmnet_all,training_data,offset_full)
  names(model_fits) <- c('HB_model', 'glmnet_model', 'training_data', 'offset')
  return(model_fits)
}