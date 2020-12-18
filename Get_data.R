Get_data <- function(wd, species, thinning, target_n_obs) {
  
  #read data if a data file has been created earlier
  thin_mark <- ifelse(thinning, "thin", "not_thin")
  if(file.exists(paste(wd, "/Compiled_data/", species, '_', thin_mark, "_quad_n_",
                       target_n_obs, '.R', sep = ''))) {
    
    load(paste(wd, "/Compiled_data/", species, '_', thin_mark, "_quad_n_",
               target_n_obs, '.R', sep = ''))
  } else {
    
    #get environmental data set
    source(paste(wd, "R_code/Get_env_data.R", sep = '/'))
    env_data <- Get_env_data(wd, species)
    
    #get species observations
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
    
    #get species data from the inventories
    source(paste(.wd, "R_code/Get_validation_data.R", sep = '/'))
    inventory_data <- Get_validation_data(wd, training_data, offset_full, species)
    
    #get environmental data for predicting
    source(paste(.wd, "R_code/Get_prediction_data.R", sep = '/'))
    pred_data <- Get_prediction_data(wd, training_data, species, scale_out)
    print("Pred_data derived")
    
    #get validation data for predicting
    source(paste(.wd, "R_code/Get_validation_data.R", sep = '/'))
    pred_PA_data <- Get_validation_data(wd, training_data, offset_full, species)
    print("Pred_PA_data derived")
    
    #derive the projection of the rasters (needed in setting folds in cross-validation)
    raster_cov_temp <- stack(paste(wd,"Data/Environment/Chelsa_SA.tif", sep = '/'))
    proj_raster <- projection(raster_cov_temp)
    
    data <- list(env_data, species_data, training_data, offset_full, inventory_data,
                 pred_data, pred_PA_data, proj_raster, scale_out)
    names(data) <- c('env_data', 'species_data', 'training_data', 'offset_full', 'inventory_data',
                     'pred_data', 'pred_PA_data', 'proj_raster', 'scale_out')
    
    save(data, file = paste(wd, "/Compiled_data/", species, '_', thin_mark, "_quad_n_",
                     target_n_obs, '.R', sep = ''))    
  }
  return(data)
}
