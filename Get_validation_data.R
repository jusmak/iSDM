Get_validation_data <- function(wd, training_data, species) {
  #Load independent validation data
  #Survey locations, get coordinates of the survey sites
  survey_sites <- read.csv(paste(wd, '/Data/PA_observations/Juan_parra_checklists/Sites_8Feb2011.csv', sep = ''))
  coord_survey_sites <- cbind(survey_sites$LongDecDeg, survey_sites$LatDecDeg)
  coord_survey_sites <- SpatialPoints(coord_survey_sites, proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  #Transform coordinates of the survey sites, get coordinate system from the raster
  raster_cov_temp <- stack(paste(wd,"Data/Environment/Chelsa_Americas.tif", sep = '/'))
  proj_temp <- projection(raster_cov_temp)
  coord_survey_sites <- spTransform(coord_survey_sites, proj_temp)
  coord_survey_sites <- as.data.frame(coord_survey_sites)
  
  #Get covariate values in the survey sites
  cov_pred <- as.data.frame(raster::extract(raster_cov_temp,coord_survey_sites))
  
  #Standardize covariate values
  cov_pred <- mapply(function(x,y,z) (x-y)/z, cov_pred, training_data$cov_mean, training_data$cov_sd)
  
  #Define locations where species is present
  survey_species <- read.csv(paste(wd, '/Data/PA_observations/Juan_parra_checklists/SpeciesxSite8Feb2011.csv', sep = ''))
  survey_species <- survey_species[,c(1:2,5)]
  #Join species with the survey locations
  survey_jp <- left_join(survey_sites, survey_species)
  
  #Define locations where species is present
  presence_loc <- data.frame(survey_jp[survey_jp$Spname=='Ocreatus_underwoodii',2])
  ind_pres <- apply(presence_loc,1, function(x) which(survey_sites[,2] == x))
  
  if (length(ind_pres == 0)) {
    
    ## Put here an empty list ##
    
  } else {
    
    #define a response variable
    resp <- rep(0, nrow(cov_pred))
    resp[ind_pres] <- 1
    
    ##WEIGHTS##
    ##EXPERT RANGE MAP##
    #load the expert prior
    expert_prior <- raster(paste(wd, "/Data/Range_map/", species, ".tif", sep = ''))
    
    #extract expert prior value in validation locations
    training_expert_prior <- raster::extract(expert_prior, coord_survey_sites)
    
    #define the lower and upper bounds of expert information
    offset_parameters <- read.csv(paste(wd, '/Data/Output_from_Diegos_draft/SDM_offset_parameters.csv', sep = ''))
    p_in <- offset_parameters[offset_parameters$ï..Species==species,2]
    #decay rate / skew / shift
    r <- offset_parameters[offset_parameters$ï..Species==species,3]
    skew <- offset_parameters[offset_parameters$ï..Species==species,4]
    shift <- offset_parameters[offset_parameters$ï..Species==species,5]
    
    #define weights
    training_expert_prior[training_expert_prior==1] <- p_in
    training_expert_prior[is.na(training_expert_prior)] <- 1-p_in
    
    #extract expert prior value in training locations
    #transform coordinates first to their original scale
    expert_range <- raster::extract(expert_prior, training_data$coordinates)
    #take only points which are inside of the range
    expert_range_inside <- training_data$coordinates[!is.na(expert_range),]
    
    #define distance from the outside validation points to the expert range
    ind_NA <- matrix(which(is.na(training_expert_prior)))
    
    
    
    
    #smooth expert prior with a convolution kernel
    #compute distance from validation points to the range edge
    dist_range <- rep(NA,nrow(ind_NA))
    for (i in 1:nrow(ind_NA)) {
      dist_range[i] <- min(sqrt(apply((training_data$coordinates[ind_NA[i],]-training_data$coordinates[-ind_NA,])^2,1,sum)))
    }
    
    #compute a smoothed prior value
    smooth_prior <- u-(u-l)/(1+exp(-r*(dist_range-shift)))^(1/skew)
    
    #combine with all other prior values
    offset_expert[is.na(training_expert_prior)] <- smooth_prior
    
    ##ELEVATION##
    #read elevation offset
    elev_range <- read.csv(paste(wd, "/Data/Elevation_range/expert_elevation.csv", sep = ''))
    elev_range_sp <- as.matrix(elev_range[elev_range$species == species,1:2])
    
    #elevation raster
    elev <- raster(paste(wd, "/Data/Elevation_range/elevation.tif", sep = ''))
    
    #elevation in the quadrature points
    elev_temp <- raster::extract(elev, xy_matrix)
    
    #indexes of quadrature points outside of the range
    ind_out <- which(elev_temp < elev_range_sp[1] | elev_temp > elev_range_sp[2])
    
    dist_elev <- matrix(rep(NA, length(ind_out)*2), nrow = length(ind_out))
    for (i in 1:length(ind_out)) {
      dist_elev[i,] <- sqrt((c(elev_temp[ind_out[i]]-elev_range_sp[1], elev_temp[ind_out[i]]-elev_range_sp[2]))^2)
      #dist_elev[i,] <- c(elev_temp[ind_out[i]]-elev_range_sp[1], elev_temp[ind_out[i]]-elev_range_sp[2])
    }
    
    #the parameter values are drawn from Diego's ms
    #decay rate / skew / shift / upper and lower asymptotes
    r <- .89
    skew <- 1
    shift <- 0
    
    #index for the cells below / above the elevation range
    ind_below <- which(elev_temp[ind_out] < elev_range_sp[1])
    ind_above <- which(elev_temp[ind_out] > elev_range_sp[2])
    
    smooth_elev_prior <- rep(NA, nrow(dist_elev))
    u <- 1
    l <- .01
    #compute a smoothed prior value for points below
    smooth_elev_prior[ind_below] <- u-(u-l)/(1+exp(-r*(dist_elev[ind_below,1]-shift)))^(1/skew)
    
    #compute a smoothed prior value for points above
    smooth_elev_prior[ind_above] <- u-(u-l)/(1+exp(-r*(dist_elev[ind_above,2]-shift)))^(1/skew)
    
    #combine elevation priors
    offset_elevation <- rep(1, length(elev_temp))
    offset_elevation[ind_out] <- smooth_elev_prior
    
    #bind all different offset scenarios
    offset_pred <- list(matrix(rep(1,length(offset_expert)*2),nrow = length(offset_expert)),
                        matrix(c(offset_expert, rep(1,length(offset_expert))),nrow = length(offset_expert)),
                        matrix(c(rep(1,length(offset_expert)), offset_elevation),nrow = length(offset_expert)),
                        matrix(c(offset_expert, offset_elevation),nrow = length(offset_expert)))
  }
  
  }
  
  
  
  
  
  
  
  
  
}