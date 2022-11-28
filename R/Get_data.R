Get_data = function(wd, species, thinning, target_n_obs, weights_area, get_offsets) {
  
  #read data if a data file has been created earlier
  thin_mark = ifelse(thinning, "thin", "not_thin")
  weights_mark = ifelse(weights_area, "wa", "not_wa")
  if (file.exists(paste0('Compiled_data/', species, '_', thin_mark, '_quad_n_',
                       target_n_obs, '_', weights_mark, '.R', sep = ''))) {
    
    load(paste0('Compiled_data/', species, '_', thin_mark, '_quad_n_',
               target_n_obs, '_', weights_mark, '.R', sep = ''))
  } else {
    
    
    #set domain
    source('R/Get_domain.R')
    Get_domain(wd, species)
    
    #get environmental data set
    source('R/Get_env_data.R')
    env_data = Get_env_data(wd, species, get_offsets)
    
    #get species po observations
    source('R/Get_species_po_data.R')
    species_po_data = Get_species_po_data(species,env_data$min_coordinates,thinning)
    
    #set scale out according to the target_n_obs
    domain_temp = stack(paste0('Data/Domains_new/', species, '.tif'))
    scale_out = round(sqrt(sum(!is.na(values(domain_temp)))/target_n_obs))
    
    #sort data for inference
    source('R/Get_species_po_training_data.R')
    PO_training_data = Get_species_po_training_data(species, env_data, species_po_data, scale_out, weights_area)
    
    #get species data from the inventories
    source('R/Get_species_pa_training_data.R')
    PA_training_data = Get_species_pa_training_data(wd, env_data, PO_training_data, species)
    
    #derive the projection of the rasters (needed in setting folds in cross-validation)
    raster_cov_temp = stack('Data/Environment/Chelsa_SA.tif')
    proj_raster = projection(raster_cov_temp)
    
    data = list(env_data, PO_training_data, PA_training_data, proj_raster, scale_out)
    names(data) = c('env_data', 'PO_training_data', 'PA_training_data', 'proj_raster', 'scale_out')
    
    save(data, file = paste0('Compiled_data/', species, '_', thin_mark, '_quad_n_',
                            target_n_obs, '_', weights_mark, '.R'))
  }
  return(data)
}
