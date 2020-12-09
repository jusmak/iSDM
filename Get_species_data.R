Get_species_data <- function(wd, species, min_coordinates, thinning) {

  #thinning and cleaning function
  source(paste(wd, "R_code/occurrence.R", sep = '/'))
  
  #define input data for thinning
  domain_temp <- stack(paste(wd,"/Data/Domains/", species, ".tif", sep = ''))
  occ_folder <- paste(wd, "/Data/PO_observations/", species, ".csv", sep = '')
  sp_occ_thin <- cleanOcc(speciesCSV = occ_folder, env = domain_temp, doThin = thinning)
  
  #derive points
  sp_occ <- as.data.frame(sp_occ_thin)[,3:4]
  
  #scale the coordinates of the species occurrence records with the same coordinates
  sp_occ <- cbind(sp_occ[,1]-min_coordinates[1],sp_occ[,2]-min_coordinates[2])
  sp_occ <- sp_occ/1000
  colnames(sp_occ) <- c('X', 'Y')

  return(sp_occ)
}