Get_species_data <- function(wd, species, min_coord) {

  # species occurrence
  sp_occ <- read.csv(paste(wd, "/Data/PO_observations/", species, ".csv", sep = ''))
  names(sp_occ) <- c('X', 'Y')
  
  #scale the species coordinates along
  sp_occ <- cbind(sp_occ[,1]-min_coord[1],sp_occ[,2]-min_coord[2])
  sp_occ <- sp_occ/1000
  
  return(sp_occ)
}