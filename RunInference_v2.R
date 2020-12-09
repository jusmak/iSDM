RunInference_v2 <- function(wd, species, observations, geographic_extent,thinning, scale_out) {
  
  # get environmental data set
  source(paste(wd, "R_code/Get_env_data.R", sep = '/'))
  env_data <- Get_env_data(wd, geographic_extent)
  
  # get species observations
  source(paste(wd, "R_code/Get_species_data.R", sep = '/'))
  species_data <- Get_species_data(wd,species,env_data$min_coordinates,thinning)
  
  #sort data for inference
  source(paste(wd, "R_code/PP_training_data.R", sep = '/'))
  training_data <- PP_training_data(wd, env_data, species_data, scale_out)
  
  #define offsets (expert range map + elevation)
  source(paste(wd, "R_code/Get_offsets.R", sep = '/'))
  offset_full <- Get_offsets(wd, training_data, env_data, species)

  #conduct Hierarchical Bayesian inference for each different offset scenario
  source(paste(wd, "R_code/Inference_HB.R", sep = '/'))
  fit_hb_all <- list()
#  for (i in 1:4) {
#    fit_hb <- Inference_HB(wd, training_data, offset_full[[i]])
#    fit_hb_all[[i]] <- fit_hb
#  }
  
  #conduct frequentist regularized regression
  source(paste(wd, "R_code/Inference_glm.R", sep = '/'))
  fit_glm_all <- list()
  for (i in 1:4) {
    fit_glm <- Inference_glm(wd, training_data, offset_full[[i]])
    fit_glm_all[[i]] <- fit_glm
  }
  
  model_fits <- list(fit_hb_all,fit_glm_all,training_data,offset_full)
  names(model_fits) <- c('HB_model', 'glm_model', 'training_data', 'offset')
  return(model_fits)
}