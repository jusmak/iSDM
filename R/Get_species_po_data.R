Get_species_po_data = function(species, min_coordinates, thinning) {

  #thinning and cleaning function
  source('R/occurrence.R')
  
  #get domain
  domain_temp = stack(paste0('Data/Domains_new/', species, '.tif'))
  
  #transform the coordinates of the records to plane coordinate system
  new_table = read.csv(paste0('Data/PO_observations_non_thinned/', species, '.csv'))
  colnames(new_table) = c('lat', 'lon')
  new_table = new_table[,c(2,1)]
  
  coordinates(new_table) = c('lon', 'lat')
  proj4string(new_table) = CRS('+init=epsg:4326') # WGS 84
  CRS.new = CRS(projection(domain_temp))
  po_coordinates_sp = spTransform(new_table, CRS.new)
  po_coordinates = po_coordinates_sp@coords
  
  #save the coordinates to a csv
  write.csv(po_coordinates, file = paste0('Data/PO_observations_non_thinned/PO_transformed_', species, '.csv'))
  
  #define input data for thinning
  occ_folder = paste0('Data/PO_observations_non_thinned/PO_transformed_', species, '.csv')
  sp_occ_thin = cleanOcc(speciesCSV = occ_folder, env = domain_temp, doThin = thinning)
  
  #derive points
  if (thinning) {
    col_ind = match(c('pres.Longitude', 'pres.Latitude'),colnames(as.data.frame(sp_occ_thin)))
  } else {
    col_ind = match(c('pres.lon', 'pres.lat'),colnames(as.data.frame(sp_occ_thin)))
  }
  
  
  sp_occ = as.data.frame(sp_occ_thin)[,col_ind]
  
  #scale the coordinates of the species occurrence records with the same coordinates
  sp_occ = cbind(sp_occ[,1]-min_coordinates[1],sp_occ[,2]-min_coordinates[2])
  sp_occ = sp_occ/1000
  colnames(sp_occ) = c('X', 'Y')

  return(sp_occ)
}