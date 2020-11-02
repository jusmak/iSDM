RunInference <- function(wd, species, offset, observations, geographic_extent) {
  
  # get environmental data set
  source(paste(wd, "R_code/Get_env_data.R", sep = '/'))
  env_data <- Get_env_data(wd, geographic_extent)
  
  # get species observations
  source(paste(wd, "R_code/Get_species_data.R", sep = '/'))
  species_data <- Get_species_data(wd,species,env_data$min_coordinates)
  
  #sort data for inference and conduct inference
  source(paste(wd, "R_code/Inference_HB.R", sep = '/'))
  
  
}