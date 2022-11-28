Get_species_po_training_data = function( species, env_data, species_po_data, scale_out, weights_area) {
  
  #define covariate values for presence-observations
  #locate them to the closest cell coordinates
  domain_temp = stack(paste0('Data/Domains_new/', species, '.tif'))
  raster_res = res(domain_temp)/1000
  
  #create a grid for covariate data
  grid_sparse = meshgrid(seq(raster_res[1]*scale_out/2, max(env_data$coordinates[,1]) - raster_res[1]*scale_out/2, raster_res[1]*scale_out),
                          seq(raster_res[2]*scale_out/2, max(env_data$coordinates[,2]) - raster_res[2]*scale_out/2, raster_res[2]*scale_out))
  
  #transform the grid into two vectors of coordinates
  quad_coarse_coordinates = cbind(array(grid_sparse$X), array(grid_sparse$Y))
  
  #collect all coordinates
  train_coordinates = rbind(species_po_data,quad_coarse_coordinates)
  
  #scale coordinates
  train_coordinates_temp = cbind(train_coordinates[,1]*1000 + env_data$min_coordinates[1],
                                  train_coordinates[,2]*1000 + env_data$min_coordinates[2])
  
  #assign covariate values to coarse grid cells
  #crop rasters with domain
  #load environmental covariates as stacked rasters
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
  
  #take covariate value only from cells which have a defined value
  values(raster_cov_temp)[-env_data$ind_na,] = NA
  #get covariate values in the presence and baclground (mesh) locations
  train_covariates = raster::extract(raster_cov_temp, train_coordinates_temp)  
  
  #find rows where any covariate has NA value
  ind_na_2 = apply(train_covariates,1,anyNA)
  train_coordinates = train_coordinates[ind_na_2==F,]
  train_covariates = train_covariates[ind_na_2==F,]
  
  # add second order effects in the training covariates
  train_covariates = cbind(train_covariates,train_covariates^2)
  colnames(train_covariates) = c(colnames(env_data$covariates), paste0(colnames(env_data$covariates), '_sqrd'))
  
  #standardize the covariates
  cov_mean = apply(train_covariates,2,mean)
  cov_sd = apply(train_covariates,2,sd)
  train_covariates_st = apply(train_covariates, 2, function(x) (x-mean(x))/sd(x))
  
  #create a response variable
  response =c(rep(1,nrow(species_po_data)),rep(0,nrow(quad_coarse_coordinates)))
  
  #assign weights to the quadarature weights and the presence observations
  presence_weights = rep(0,nrow(species_po_data))
  
  #weights of the fine scale quadrature points
  #quad_fine_weights = rep(1, nrow(quad_fine_coordinates))
  
  #weights for sparse quadrature points are scale_out^2 or 1
  weight_coarse = ifelse(weights_area, scale_out^2, 1)
  quad_coarse_weights = rep(weight_coarse, nrow(quad_coarse_coordinates))
  
  #combine all weights into one vector
  weights = c(presence_weights, quad_coarse_weights)
  
  #subsample response and weights
  weights = weights[ind_na_2==F]
  response = response[ind_na_2==F]

  #collect offset information in the cells
  #temp coordinates for picking offset values from raster
  temp_coord = cbind(train_coordinates[,1]*1000 + env_data$min_coordinates[1], train_coordinates[,2]*1000 + env_data$min_coordinates[2])

  #create a raster template with all cells NA
  temp_raster = domain_temp
  values(temp_raster) = NA
  
  #expert map
  #define values in raster
  values(temp_raster)[env_data$ind_na] = env_data$offsets[,1]
  #extract values in data points
  offset_expert = raster::extract(temp_raster, temp_coord)
  offset_expert[is.na(offset_expert)] = min(offset_expert, na.rm = TRUE)
  #scale offsets to sum to the area of the study points
  pop_i = sum(offset_expert*weights)
  offset_expert = offset_expert/pop_i
  #center the offset to zero
  offset_expert = offset_expert*sum(weights)
  
  #elevation limits
  #define values in raster
  values(temp_raster)[env_data$ind_na] = env_data$offsets[,2]
  #extract values in quadrature points
  offset_elevation = raster::extract(temp_raster, temp_coord)
  offset_elevation[is.na(offset_elevation)] = min(offset_elevation, na.rm = TRUE)
  #scale offsets to sum to the area of the study points
  pop_i = sum(offset_expert*weights)
  offset_elevation = offset_elevation/pop_i
  #center the offset to zero
  offset_elevation = offset_elevation*sum(weights)
  
  #define training data as a list
  training_data = list(train_coordinates, train_covariates_st, cov_mean, cov_sd, env_data$min_coordinates, response, weights, offset_expert, offset_elevation)
  names(training_data) = c('coordinates', 'covariates', 'cov_mean', 'cov_sd', 'min_coordinates', 'response', 'weights', 'offset_expert', 'offset_elevation')
  
  return(training_data)
}