Get_env_data = function(wd, species, get_offsets){
  
  #load environmental covariates as stacked rasters
  raster_cov_temp = stack('Data/Environment/Chelsa_SA.tif')
  elevation_temp = stack('Data/Elevation_range/topography_elevation_1KMmi_GMTEDmi_NA30x30_americas.tif')
  cloud_temp = raster('Data/Cloud_cover/MODCF_intraannualSD_resampled_masked_NA30x30_americas.tif')
  evi_temp = raster('Data/EVI/Annual_EVI_resampled_NA30x30_americas.tif')
  tri_temp = raster('Data/TRI/tri_1KMmd_GMTEDmd_resampled_masked_NA30x30_americas.tif')
  
  #load geographic extent
  domain_temp = stack(paste0('Data/Domains_new/', species, '.tif'))
  
  #combine all raster layers into one stack
  raster_cov_temp = stack(raster_cov_temp, cloud_temp, evi_temp, tri_temp)  
  
  #crop covariates with the geographic extent
  raster_cov_temp = crop(raster_cov_temp, domain_temp)
  
  #turn all cells not defined in domain to NA
  values(raster_cov_temp)[is.na(values(domain_temp)),] = NA
  
  #rescale covariates to 25x25km cell size
  raster_cov_temp = aggregate(raster_cov_temp, fact = 25)
  
  #find cells which have all covariate values
  cov_df = as.matrix(raster_cov_temp, na.rm = FALSE)
  rows_na = apply(cov_df,1,anyNA)
  ind_na = which(rows_na==FALSE)
  cov_matrix = cov_df[ind_na,]
  
  #get coordinates of the cells
  xy_matrix = coordinates(raster_cov_temp)[ind_na,]
  
  ## Get offsets
  offsets = matrix(1, nrow = nrow(xy_matrix), ncol = 2)
  
  if (get_offsets == T) {
  ##Here spatial offsets##
  }

  ## Collect all data
  #put coordinates start from zero
  min_coord = apply(xy_matrix,2,min)
  xy_matrix = cbind(xy_matrix[,1]-min_coord[1],xy_matrix[,2]-min_coord[2])
  
  #scale coordinates to km:s
  xy_matrix = xy_matrix/1000
  
  #assign correct names for both matrices
  colnames(cov_matrix) = c('cov1', 'cov2', 'cov3', 'cov4', 'cloud', 'evi', 'tri')
  colnames(xy_matrix) = c('X', 'Y')
  names(min_coord) = c('X', 'Y')

  #list matrices
  all_data = list(cov_matrix,xy_matrix,min_coord,offsets,ind_na)
  names(all_data) = c('covariates', 'coordinates', 'min_coordinates','offsets','ind_na')
  
  return(all_data)
}