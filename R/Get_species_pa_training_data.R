Get_species_pa_training_data = function( env_data, PO_training_data, species) {
  
  #Load independent validation data
  #Survey locations, get coordinates of the survey sites
  survey_sites = read.csv('Data/PA_observations/Juan_parra_checklists/Sites_8Feb2011.csv')
  coord_survey_sites = cbind(survey_sites$LongDecDeg, survey_sites$LatDecDeg)
  coord_survey_sites = SpatialPoints(coord_survey_sites, proj4string=CRS('+proj=longlat +datum=WGS84'))
  
  #Transform coordinates of the survey sites, get coordinate system from the raster
  raster_cov_temp = stack('Data/Environment/Chelsa_SA.tif')
  proj_temp = projection(raster_cov_temp)
  coord_survey_sites = spTransform(coord_survey_sites, proj_temp)
  coord_survey_sites = as.data.frame(coord_survey_sites)
  
  #Take only the survey sites which are inside of the domain
  domain_temp = stack(paste0('/Data/Domains_new/', species, '.tif'))
  ind_site = raster::extract(domain_temp, coord_survey_sites)
  coord_survey_sites = coord_survey_sites[!is.na(ind_site),]
  
  #load environmental covariates as stacked rasters
  raster_cov_temp = stack('Data/Environment/Chelsa_SA.tif')
  elevation_temp = stack('Data/Elevation_range/topography_elevation_1KMmi_GMTEDmi_NA30x30_americas.tif')
  cloud_temp = raster('Data/Cloud_cover/MODCF_intraannualSD_resampled_masked_NA30x30_americas.tif')
  evi_temp = raster('Data/EVI/Annual_EVI_resampled_NA30x30_americas.tif')
  tri_temp = raster('Data/TRI/tri_1KMmd_GMTEDmd_resampled_masked_NA30x30_americas.tif')
  
  #load geographic extent
  domain_temp = stack(paste0('/Data/Domains_new/', species, '.tif'))
  
  #combine all raster layers into one stack
  raster_cov_temp = stack(raster_cov_temp, cloud_temp, evi_temp, tri_temp)  
  
  #crop covariates with the geographic extent
  raster_cov_temp = crop(raster_cov_temp, domain_temp)
  
  #Get covariate values in the survey sites
  cov_temp = as.data.frame(raster::extract(raster_cov_temp,coord_survey_sites))
  
  #Get second order effects
  train_covariates = cbind(cov_temp, cov_temp^2)
  
  #Standardize covariate values
  train_covariates_st = matrix(NA, nrow = nrow(train_covariates), ncol = ncol(train_covariates))
  for (j in 1:ncol(train_covariates)) {
    train_covariates_st[,j] = (train_covariates[,j]-PO_training_data$cov_mean[j])/PO_training_data$cov_sd[j]
  }
  
  #train_covariates_st = mapply(function(x,y,z) (x-y)/z, train_covariates, PO_training_data$cov_mean, PO_training_data$cov_sd)
  
  #Define column names
  colnames(train_covariates_st) = c(colnames(env_data$covariates), paste0(colnames(env_data$covariates), '_sqrd'))
  
  #Define locations where species is present
  survey_species = read.csv('Data/PA_observations/Juan_parra_checklists/SpeciesxSite8Feb2011.csv')
  survey_species = survey_species[,c(1:2,5)]
  pres_sp = survey_species[survey_species$Spname==species,]
  
  #Join presences with survey locations
  survey_jp = left_join(survey_sites[!is.na(ind_site),], pres_sp)
  
  #define a response variable
  resp = rep(0, nrow(train_covariates_st))
  resp[!is.na(survey_jp$Spname)] = 1
  
  #transform coordinates
  coord_pred = cbind(coord_survey_sites[,1]-PO_training_data$min_coordinates[1],
                      coord_survey_sites[,2]-PO_training_data$min_coordinates[2])/1000
  
  colnames(coord_pred) = c('X', 'Y')

  ##EXPERT RANGE MAP and ELEVATION##
  #get exact values from the raster
  #create a raster template with all cells NA
  temp_raster = domain_temp
  values(temp_raster) = NA
  
  #expert map
  #define values in raster
  values(temp_raster)[env_data$ind_na] = env_data$offsets[,1]
  #extract values in data points
  offset_expert = raster::extract(temp_raster, coord_survey_sites)
  offset_expert[is.na(offset_expert)] = min(offset_expert, na.rm = TRUE)
  #scale offsets according to the presence-only points
  pop_i = sum(PO_training_data$offset_expert*PO_training_data$weights)
  offset_expert = offset_expert/pop_i
  #center the offset to zero
  offset_expert = offset_expert*sum(PO_training_data$weights)
  
  #elevation
  #define values in raster
  values(temp_raster)[env_data$ind_na] = env_data$offsets[,2]
  #extract values in data points
  offset_elevation = raster::extract(temp_raster, coord_survey_sites)
  offset_elevation[is.na(offset_elevation)] = min(offset_elevation, na.rm = TRUE)
  #scale offsets according to the presence-only points
  pop_i = sum(PO_training_data$offset_elevation*PO_training_data$weights)
  offset_elevation = offset_elevation/pop_i
  #center the offset to zero
  offset_elevation = offset_elevation*sum(PO_training_data$weights)
  
  
  survey_data = list(coord_pred, train_covariates_st, resp, offset_expert, offset_elevation)
  names(survey_data) = c('coordinates', 'covariates', 'response', 'offset_expert', 'offset_elevation')
  
  return(survey_data)
}