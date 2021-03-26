Get_validation_data <- function(wd, training_data, offset, species) {
  #Load independent validation data
  #Survey locations, get coordinates of the survey sites
  survey_sites <- read.csv(paste(wd, '/Data/PA_observations/Juan_parra_checklists/Sites_8Feb2011.csv', sep = ''))
  coord_survey_sites <- cbind(survey_sites$LongDecDeg, survey_sites$LatDecDeg)
  coord_survey_sites <- SpatialPoints(coord_survey_sites, proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  #Transform coordinates of the survey sites, get coordinate system from the raster
  raster_cov_temp <- stack(paste(wd,"Data/Environment/Chelsa_SA.tif", sep = '/'))
  proj_temp <- projection(raster_cov_temp)
  coord_survey_sites <- spTransform(coord_survey_sites, proj_temp)
  coord_survey_sites <- as.data.frame(coord_survey_sites)
  
  #Take only the survey sites which are inside of the domain
  domain_temp <- stack(paste(wd,"/Data/Domains/", species, ".tif", sep = ''))
  ind_site <- raster::extract(domain_temp[[2]], coord_survey_sites)
  coord_survey_sites <- coord_survey_sites[!is.na(ind_site),]
  
  #Get covariate values in the survey sites
  cov_pred <- as.data.frame(raster::extract(raster_cov_temp,coord_survey_sites))
  
  #Standardize covariate values
  cov_pred <- mapply(function(x,y,z) (x-y)/z, cov_pred, training_data$cov_mean, training_data$cov_sd)
  
  #Define locations where species is present
  survey_species <- read.csv(paste(wd, '/Data/PA_observations/Juan_parra_checklists/SpeciesxSite8Feb2011.csv', sep = ''))
  survey_species <- survey_species[,c(1:2,5)]
  pres_sp <- survey_species[survey_species$Spname==species,]
  
  #Join presences with survey locations
  survey_jp <- left_join(survey_sites[!is.na(ind_site),], pres_sp)
  
  #define a response variable
  resp <- rep(0, nrow(cov_pred))
  resp[!is.na(survey_jp$Spname)] <- 1
  
  #transform coordinates
  coord_pred <- cbind(coord_survey_sites[,1]-training_data$min_coordinates[1],
                      coord_survey_sites[,2]-training_data$min_coordinates[2])/1000

  ##EXPERT RANGE MAP and ELEVATION##
  #find the closest training point for each validation point
  ind_tr <- rep(NA,length(resp))
  for (i in 1:length(ind_tr)) {
    ind_tr[i] <- which.min((coord_pred[i,1]-training_data$coordinates[,1])^2 + 
                             (coord_pred[i,2]-training_data$coordinates[,2])^2)
  }
  
  #define offsets for each validation point
  offset_pred <- list()
  for (i in 1:length(offset)) {
    offset_pred[[i]] <- offset[[i]][ind_tr,]
  }
  
  independent_validation <- list(coord_pred, cov_pred, resp, offset_pred)
  names(independent_validation) <- c('coordinates', 'covariates', 'response', 'offset')
  
  return(independent_validation)
}