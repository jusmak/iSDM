Get_data <- function(wd, species, thinning, target_n_obs, weights_area) {
  
  #read data if a data file has been created earlier
  thin_mark <- ifelse(thinning, "thin", "not_thin")
  weights_mark <- ifelse(weights_area, "wa", "not_wa")
  if(file.exists(paste(wd, "/Compiled_data/", species, '_', thin_mark, "_quad_n_",
                       target_n_obs, "_", weights_mark, '.R', sep = ''))) {
    
    load(paste(wd, "/Compiled_data/", species, '_', thin_mark, "_quad_n_",
               target_n_obs, "_", weights_mark, '.R', sep = ''))
  } else {
    
    #get environmental data set
    source(paste(wd, "R_code/Get_env_data.R", sep = '/'))
    env_data <- Get_env_data(wd, species)
    
    #get species po observations
    source(paste(wd, "R_code/Get_species_po_data.R", sep = '/'))
    species_po_data <- Get_species_po_data(wd,species,env_data$min_coordinates,thinning)
    
    #set scale out according to the target_n_obs
    domain_temp <- stack(paste(wd,"/Data/Domains/", species, ".tif", sep = ''))
    scale_out <- round(sqrt(sum(!is.na(values(domain_temp[[2]])))/target_n_obs))
    
    #sort data for inference
    source(paste(wd, "R_code/PP_training_data.R", sep = '/'))
    training_data <- PP_training_data(wd, species, env_data, species_po_data, scale_out, weights_area)
    
    #get species data from the inventories
    source(paste(.wd, "R_code/Get_species_pa_data.R", sep = '/'))
    inventory_data <- Get_species_pa_data(wd, env_data, training_data, species)
    
    #derive the projection of the rasters (needed in setting folds in cross-validation)
    raster_cov_temp <- stack(paste(wd,"Data/Environment/Chelsa_SA.tif", sep = '/'))
    proj_raster <- projection(raster_cov_temp)
    
    data <- list(env_data, training_data, inventory_data, proj_raster, scale_out)
    names(data) <- c('env_data', 'training_data', 'inventory_data',
                     'proj_raster', 'scale_out')
    
    save(data, file = paste(wd, "/Compiled_data/", species, '_', thin_mark, "_quad_n_",
                            target_n_obs, "_", weights_mark, '.R', sep = ''))    
  }
  return(data)
}
