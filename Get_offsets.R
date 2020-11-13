Get_offsets <- function(wd, training_data, env_data, species) {
  
  ##EXPERT RANGE MAP##
  #load the expert prior
  expert_prior <- raster(paste(wd, "/Data/Range_map/", species, ".tif", sep = ''))
  
  #extract expert prior value in training locations
  coord_temp <- training_data$coordinates*1000
  coord_temp <- cbind(coord_temp[,1]+env_data$min_coordinates[1], coord_temp[,2]+env_data$min_coordinates[2])
  training_expert_prior <- raster::extract(expert_prior, coord_temp)
  
  #define the lower and upper bounds of expert information
  offset_parameters <- read.csv(paste(wd, '/Data/Output_from_Diegos_draft/SDM_offset_parameters.csv', sep = ''))
  p_in <- offset_parameters[offset_parameters$ï..Species==species,2]
  #decay rate / skew / shift
  r <- offset_parameters[offset_parameters$ï..Species==species,3]
  skew <- offset_parameters[offset_parameters$ï..Species==species,4]
  shift <- offset_parameters[offset_parameters$ï..Species==species,5]
  
  offset_expert <- training_expert_prior
  #assign probability for the cells inside of the range
  offset_expert[!is.na(training_expert_prior)] <- (p_in)*sum(training_data$response)/sum(!is.na(training_expert_prior))
  #assign probability for the cells outside of the range
  offset_expert[is.na(training_expert_prior)] <- (1-p_in)*sum(training_data$response)/sum(is.na(training_expert_prior))
  #upper / lower asymptote
  u <- (p_in)/sum(!is.na(training_expert_prior))
  l <- (1-p_in)/sum(is.na(training_expert_prior))
  
  #smooth expert prior with a convolution kernel
  #compute distance to the range edge for the cells outside of the range
  ind_NA <- matrix(which(is.na(training_expert_prior)))
  
  dist_range <- rep(NA,nrow(ind_NA))
  for (i in 1:nrow(ind_NA)) {
    dist_range[i] <- min(sqrt(apply((training_data$coordinates[ind_NA[i],]-training_data$coordinates[-ind_NA,])^2,1,sum)))
  }
  
  #compute a smoothed prior value
  smooth_prior <- u-(u-l)/(1+exp(-r*(dist_range-shift)))^(1/skew)
  
  #normalize the prior values to equal the original sum of probabilities outside of the range
  smooth_prior_sc <- smooth_prior/(sum(smooth_prior)/((1-p_in)*sum(training_data$response)))
  
  #combine with all other prior values
  offset_expert[is.na(training_expert_prior)] <- smooth_prior_sc
  
  ##ELEVATION##
  #read elevation offset
  elev_range <- read.csv(paste(wd, "/Data/Elevation_range/expert_elevation.csv", sep = ''))
  elev_range_sp <- as.matrix(elev_range[elev_range$species == species,1:2])
  
  #elevation raster
  elev <- raster(paste(wd, "/Data/Elevation_range/elevation.tif", sep = ''))
  
  #elevation in the quadrature points
  elev_temp <- raster::extract(elev, coord_temp)
  
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
  u <- 100
  l <- 1
  #compute a smoothed prior value for points below
  smooth_elev_prior[ind_below] <- u-(u-l)/(1+exp(-r*(dist_elev[ind_below,1]-shift)))^(1/skew)
  
  #compute a smoothed prior value for points above
  smooth_elev_prior[ind_above] <- u-(u-l)/(1+exp(-r*(dist_elev[ind_above,2]-shift)))^(1/skew)
  
  #combine elevation priors
  offset_elevation <- rep(100, length(elev_temp))
  offset_elevation[ind_out] <- smooth_elev_prior
  
  #make it sum to the number of presence observations
  offset_elevation <- offset_elevation*sum(training_data$response)/sum(offset_elevation)
  
  #bind all different offset scenarios
  offset_full <- list(matrix(rep(sum(training_data$response)/length(offset_expert),length(offset_expert)*2),nrow = length(offset_expert)),
                      matrix(c(offset_expert, rep(sum(training_data$response)/length(offset_expert),length(offset_expert))),nrow = length(offset_expert)),
                      matrix(c(rep(sum(training_data$response)/length(offset_expert),length(offset_expert)), offset_elevation),nrow = length(offset_expert)),
                      matrix(c(offset_expert, offset_elevation),nrow = length(offset_expert)))
                      
  return(offset_full)
}
