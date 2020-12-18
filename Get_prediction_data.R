#Get prediction data
Get_prediction_data <- function(wd, training_data, species, scale_out) {
  #load raster
  raster_cov_temp <- stack(paste(wd,"Data/Environment/Chelsa_SA.tif", sep = '/'))
  
  #load geographic extent
  domain_temp <- stack(paste(wd,"/Data/Domains/", species, ".tif", sep = ''))
  raster_cov_temp <- crop(raster_cov_temp, domain_temp)
  values(raster_cov_temp)[is.na(values(domain_temp)[,2]),] <- NA
  
  #rescale raster
  raster_cov_temp <- raster::aggregate(raster_cov_temp, fact = scale_out)
  
  #name the raster layers
  names(raster_cov_temp) <- c("cov1", "cov2", "cov3", "cov4")

  #standardize covariates
  cov_df <- as.data.frame(raster_cov_temp, na.rm = FALSE)
  rows_na <- apply(cov_df,1,anyNA)
  cov_df <- cov_df[rows_na==FALSE,]
  cov_pred <- mapply(function(x,y,z) (x-y)/z, cov_df, training_data$cov_mean, training_data$cov_sd)
  
  #get coordinates of the cells
  xy_matrix <- coordinates(raster_cov_temp)
  xy_matrix <- xy_matrix[rows_na==FALSE,]
  
  #put coordinates start from zero
  xy_sc_matrix <- cbind(xy_matrix[,1]-training_data$min_coordinates[1],xy_matrix[,2]-training_data$min_coordinates[2])
  
  #scale coordinates to km:s
  xy_sc_matrix <- xy_sc_matrix/1000
  coord_pred <- xy_sc_matrix
  
  
  ##EXPERT RANGE MAP##
  #load the expert prior
  expert_prior <- raster(paste(wd, "/Data/Range_map/", species, ".tif", sep = ''))
  
  #extract expert prior value in training locations
  expert_range <- raster::extract(expert_prior, xy_matrix)
  
  #define the lower and upper bounds of expert information
  offset_parameters <- read.csv(paste(wd, '/Data/Output_from_Diegos_draft/SDM_offset_parameters.csv', sep = ''))
  colnames(offset_parameters)[1] <- 'Species'
  p_in <- offset_parameters[offset_parameters$Species==species,2]
  #decay rate / skew / shift
  r <- offset_parameters[offset_parameters$Species==species,3]
  skew <- offset_parameters[offset_parameters$Species==species,4]
  shift <- offset_parameters[offset_parameters$Species==species,5]
  
  offset_expert <- expert_range
  #assign intensity for the cells inside of the range
  offset_expert[!is.na(expert_range)] <- p_in*sum(training_data$response)/(sum(!is.na(expert_range))*scale_out^2)
  #assign intensity for the cells outside of the range
  offset_expert[is.na(expert_range)] <- (1-p_in)*sum(training_data$response)/(sum(is.na(expert_range))*scale_out^2)
  #upper / lower asymptote of the intensity
  u <- p_in*sum(training_data$response)/(sum(!is.na(expert_range))*scale_out^2)
  l <- (1-p_in)*sum(training_data$response)/(sum(is.na(expert_range))*scale_out^2)
  
  #smooth expert prior with a convolution kernel
  #compute distance to the range edge for the cells outside of the range
  #take only points which are inside the range
  expert_range_inside <- coord_pred[!is.na(expert_range),]
  
  #training points which are outside the range
  ind_NA <- matrix(which(is.na(expert_range)))
  
  dist_range <- rep(NA,nrow(ind_NA))
  for (i in 1:nrow(ind_NA)) {
    dist_range[i] <- min(sqrt((coord_pred[ind_NA[i],1]-expert_range_inside[,1])^2 + 
                                      (coord_pred[ind_NA[i],2]-expert_range_inside[,2])^2))
  }
  
  #compute a smoothed prior value
  smooth_prior <- u-(u-l)/(1+exp(-r*(dist_range-shift)))^(1/skew)

  #derive the predicted population
  pop_i <- smooth_prior * scale_out^2
  
  #scale population to equal the number of presence observations
  pop_i <- pop_i/sum(pop_i)*((1-p_in)*sum(training_data$response))
  
  #derive the intensity for each cell
  smooth_prior_sc <- pop_i/scale_out^2
  
  #combine with all other prior values
  offset_expert[is.na(expert_range)] <- smooth_prior_sc
  
  #center the offset to zero
  offset_expert <- offset_expert/(sum(training_data$response)/(scale_out^2*length(offset_expert)))
  
  ##ELEVATION##
  #read elevation offset
  elev_range <- read.csv(paste(wd, "/Data/Elevation_range/expert_elevation.csv", sep = ''))
  elev_range_sp <- as.matrix(elev_range[elev_range$species == species,1:2])
  
  #elevation raster
  elev <- raster(paste(wd, "/Data/Elevation_range/elevation.tif", sep = ''))
  
  #elevation in the prediction points
  elev_temp <- raster::extract(elev, xy_matrix)
  
  #indexes of quadrature points outside of the range
  ind_out <- which(elev_temp < elev_range_sp[1] | elev_temp > elev_range_sp[2])
  
  dist_elev <- matrix(rep(NA, length(ind_out)*2), nrow = length(ind_out))
  for (i in 1:length(ind_out)) {
    dist_elev[i,] <- sqrt((c(elev_temp[ind_out[i]]-elev_range_sp[1], elev_temp[ind_out[i]]-elev_range_sp[2]))^2)
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
  #intensity between elevation limits
  u <- .99*sum(training_data$response)/(length(elev_temp[-ind_out])*scale_out^2)
  
  #intensity outside of the elevation limits
  l <- .01*sum(training_data$response)/(length(elev_temp[ind_out])*scale_out^2)

  #compute a smoothed prior value for points below
  smooth_elev_prior[ind_below] <- u-(u-l)/(1+exp(-r*(dist_elev[ind_below,1]-shift)))^(1/skew)
  
  #compute a smoothed prior value for points above
  smooth_elev_prior[ind_above] <- u-(u-l)/(1+exp(-r*(dist_elev[ind_above,2]-shift)))^(1/skew)

  #derive the predicted population
  pop_i <- smooth_elev_prior * scale_out^2
  
  #scale population to equal the number of presence observations
  pop_i <- pop_i/sum(pop_i)*(.01*sum(training_data$response))
  
  #derive the intensity for each cell
  smooth_elev_prior_sc <- pop_i/scale_out^2
  
  #combine elevation priors
  offset_elevation <- rep(NA, length(elev_temp))
  offset_elevation[ind_out] <- smooth_elev_prior_sc
  offset_elevation[-ind_out] <- u
  
  #center the offset to zero
  offset_elevation <- offset_elevation/(sum(training_data$response)/(scale_out^2*length(offset_expert)))
  
  #bind all different offset scenarios
  offset_full <- list(matrix(rep(1,length(offset_expert)*2),nrow = length(offset_expert)),
                      matrix(c(offset_expert, rep(1,length(offset_expert))),nrow = length(offset_expert)),
                      matrix(c(rep(1,length(offset_expert)), offset_elevation),nrow = length(offset_expert)),
                      matrix(c(offset_expert, offset_elevation),nrow = length(offset_expert)))
  
  pred_data <- list(cov_pred, coord_pred, offset_full)
  names(pred_data) <- c('covariates', 'coordinates', 'offset')
  return(pred_data)
}