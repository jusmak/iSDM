Get_env_data <- function(wd, species){
  
  #load environmental covariates as stacked rasters
  raster_cov_temp <- stack(paste(wd,"Data/Environment/Chelsa_SA.tif", sep = '/'))

  #load geographic extent
  domain_temp <- stack(paste(wd,"/Data/Domains/", species, ".tif", sep = ''))
  
  #crop covariates with the geographic extent
  raster_cov_temp <- crop(raster_cov_temp, domain_temp)
  
  #turn all cells not defined in domain to NA
  values(raster_cov_temp)[is.na(values(domain_temp)[,2]),] <- NA
  
  #find cells which have all covariate values
  cov_df <- as.matrix(raster_cov_temp, na.rm = FALSE)
  rows_na <- apply(cov_df,1,anyNA)
  ind_na <- which(rows_na==FALSE)
  cov_matrix <- cov_df[ind_na,]
  
  #get coordinates of the cells
  xy_matrix <- coordinates(raster_cov_temp)[ind_na,]
  
  
  ## Get offsets
  
  ##EXPERT RANGE MAP##
  #load the expert prior
  expert_prior <- raster(paste(wd, "/Data/Range_map/", species, ".tif", sep = ''))
  
  #crop the expert map
  expert_prior <- crop(expert_prior, domain_temp)
  
  #compute number of presence and absence cells
  n_pres <- sum(!is.na(values(expert_prior)))
  n_abs <- length(ind_na)-n_pres
  
  #define the lower and upper bounds of expert information
  offset_parameters <- read.csv(paste(wd, '/Data/Output_from_Diegos_draft/SDM_offset_parameters.csv', sep = ''))
  colnames(offset_parameters)[1] <- 'Species'
  p_in <- offset_parameters[offset_parameters$Species==species,2]
  #decay rate / skew / shift
  r <- offset_parameters[offset_parameters$Species==species,3]
  skew <- offset_parameters[offset_parameters$Species==species,4]
  shift <- offset_parameters[offset_parameters$Species==species,5]
  
  #upper / lower asymptote of the intensity when continuously mapping over the domain
  u <- p_in/n_pres
  l <- (1-p_in)/n_abs
  
  #smooth expert prior with a convolution kernel
  #compute distance to the range edge for the cells outside of the range
  dist <- distance(expert_prior, doEdge=TRUE)
  dist_m <- values(dist)/1000
  
  #compute a smoothed prior value
  smooth_prior <- u-(u-l)/(1+exp(-r*(dist_m-shift)))^(1/skew)
  
  #add constant values to the cells inside the range
  smooth_prior[!is.na(values(expert_prior))] <- u
  
  #keep only the cells that have covariate values
  smooth_prior <- smooth_prior[ind_na]
  
  #derive the predicted population
  pop_i <- sum(smooth_prior)
  
  #scale intensities to match the population to the number of presence observations
  smooth_prior_sc <- smooth_prior/pop_i
  
  #center the offset to zero
  offset_expert <- smooth_prior_sc*length(smooth_prior_sc)
  
  ##ELEVATION##
  #read elevation offset
  elev_range <- read.csv(paste(wd, "/Data/Elevation_range/expert_elevation.csv", sep = ''))
  elev_range_sp <- as.matrix(elev_range[elev_range$species == species,1:2])
  
  #elevation raster
  elev <- raster(paste(wd, "/Data/Elevation_range/elevation.tif", sep = ''))
  
  #crop elevation with domain
  elev_temp <- crop(elev, domain_temp[[2]])
  elev_values <- values(elev_temp)[ind_na]
  
  #number of presences and absences
  n_pres <- sum((elev_values >= elev_range_sp[1] & elev_values <= elev_range_sp[2]))
  n_abs <- sum(elev_values < elev_range_sp[1] | elev_values > elev_range_sp[2])
  
  #indexes of data points outside of the range
  ind_out <- which(elev_values < elev_range_sp[1] | elev_values > elev_range_sp[2])
  
  dist_elev <- matrix(rep(NA, length(ind_out)*2), nrow = length(ind_out))
  for (i in 1:length(ind_out)) {
    dist_elev[i,] <- sqrt((c(elev_values[ind_out[i]]-elev_range_sp[1], elev_values[ind_out[i]]-elev_range_sp[2]))^2)
  }
  
  #the parameter values are drawn from Diego's ms
  #decay rate / skew / shift / upper and lower asymptotes
  r <- .89
  skew <- 1
  shift <- 0
  
  #index for the cells below / above the elevation range
  ind_below <- which(elev_values[ind_out] < elev_range_sp[1])
  ind_above <- which(elev_values[ind_out] > elev_range_sp[2])
  
  smooth_elev_prior <- rep(NA, nrow(dist_elev))
  
  #intensity between elevation limits
  u <- .99/n_pres
  
  #intensity outside of the elevation limits
  l <- .01/n_abs
  
  #compute a smoothed prior value for points below
  smooth_elev_prior[ind_below] <- u-(u-l)/(1+exp(-r*(dist_elev[ind_below,1]-shift)))^(1/skew)
  
  #compute a smoothed prior value for points above
  smooth_elev_prior[ind_above] <- u-(u-l)/(1+exp(-r*(dist_elev[ind_above,2]-shift)))^(1/skew)

  #combine elevation priors
  offset_elevation <- rep(NA, length(elev_values))
  offset_elevation[ind_out] <- smooth_elev_prior
  offset_elevation[-ind_out] <- u

  #derive the predicted population
  pop_i <- sum(offset_elevation)
  
  #scale offsets according to the sum
  offset_elevation <- offset_elevation/pop_i
  
  #center the offset to zero
  offset_elevation <- offset_elevation*length(offset_elevation)
  
  offsets <- cbind(offset_expert, offset_elevation)

  ## Collect all data
  
  #put coordinates start from zero
  min_coord <- apply(xy_matrix,2,min)
  xy_matrix <- cbind(xy_matrix[,1]-min_coord[1],xy_matrix[,2]-min_coord[2])
  
  #scale coordinates to km:s
  xy_matrix <- xy_matrix/1000
  
  #assign correct names for both matrices
  colnames(cov_matrix) <- c('cov1', 'cov2', 'cov3', 'cov4')
  colnames(xy_matrix) <- c('X', 'Y')
  names(min_coord) <- c('X', 'Y')

  #list matrices
  all_data <- list(cov_matrix,xy_matrix,min_coord,offsets,ind_na)
  names(all_data) <- c("covariates", "coordinates", "min_coordinates","offsets","ind_na")
  
  return(all_data)
}