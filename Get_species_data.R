Get_species_data <- function(wd, species, min_coordinates, thinning) {

#  #species occurrence
#  sp_occ <- read.csv(paste(wd, "/Data/PO_observations/", species, ".csv", sep = ''))
  
  #thinning and cleaning function
  source(paste(wd, "R_code/occurrence.R", sep = '/'))
  
  #define input data for thinning
  occ_folder <- paste(wd, "/Data/PO_observations/", species, ".csv", sep = '')
  env_raster <- raster_cov_temp <- stack(paste(wd,"Data/Environment/Chelsa_SA.tif", sep = '/'))
  sp_occ_thin <- cleanOcc(speciesCSV = occ_folder, env = env_raster, doThin = thinning)
  
  #derive points
  sp_occ <- as.data.frame(sp_occ_thin)[,3:4]
  
  #scale the species coordinates along
  sp_occ <- cbind(sp_occ[,1]-min_coordinates[1],sp_occ[,2]-min_coordinates[2])
  sp_occ <- sp_occ/1000
  colnames(sp_occ) <- c('X', 'Y')
  
  
  
  return(sp_occ)
}