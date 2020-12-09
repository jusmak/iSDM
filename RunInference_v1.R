RunInference_v1 <- function(wd, species, thinning, target_n_obs, inference_type, HB_inf, glmnet_inf, ppml_inf) {
  
  # get environmental data set
  source(paste(wd, "R_code/Get_env_data.R", sep = '/'))
  env_data <- Get_env_data(wd, species)
  
  # get species observations
  source(paste(wd, "R_code/Get_species_data.R", sep = '/'))
  species_data <- Get_species_data(wd,species,env_data$min_coordinates,thinning)
  
  #set scale out according to the target_n_obs
  domain_temp <- stack(paste(wd,"/Data/Domains/", species, ".tif", sep = ''))
  scale_out <- round(sqrt(sum(!is.na(values(domain_temp[[1]])))/50000))
  
  #sort data for inference
  source(paste(wd, "R_code/PP_training_data.R", sep = '/'))
  training_data <- PP_training_data(wd, env_data, species_data, scale_out)
  
  #define offsets (expert range map + elevation)
  source(paste(wd, "R_code/Get_offsets.R", sep = '/'))
  offset_full <- Get_offsets(wd, training_data, env_data, species)

  #conduct Hierarchical Bayesian inference for each different offset scenario
  fit_hb_all <- list()
  if (HB_inf) {
    source(paste(wd, "R_code/Inference_HB.R", sep = '/'))
    for (i in 1:4) {
      fit_hb <- Inference_HB(wd, training_data, offset_full[[i]])
      fit_hb_all[[i]] <- fit_hb
    }
  }
  
  #conduct frequentist regularized regression for each different offset scenario
  fit_glmnet_all <- list()
  if (glmnet_inf) {
    source(paste(wd, "R_code/Inference_glmnet.R", sep = '/'))
    for (i in 1:4) {
      fit_glmnet <- Inference_glmnet(wd, training_data, offset_full[[i]])
      fit_glmnet_all[[i]] <- fit_glmnet
    }
  }
  
  #conduct frequentist regularized regression with ppmlasso for each different offset scenario
  fit_ppml_all <- list()
  if (ppml_inf) {
    source(paste(wd, "R_code/Inference_ppml.R", sep = '/'))
    for (i in 1:4) {
      fit_ppml <- Inference_ppml(wd, training_data, offset_full[[i]])
      fit_ppml_all[[i]] <- fit_ppml
    }
  }
  
  model_fits <- list(fit_hb_all,fit_glmnet_all,fit_ppml_all,training_data,offset_full)
  names(model_fits) <- c('HB_model', 'glmnet_model', 'ppml_model', 'training_data', 'offset')
  return(model_fits)
}