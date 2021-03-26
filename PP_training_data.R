PP_training_data <- function(wd, species, env_data, species_data, scale_out, weights_area) {
  
  #define covariate values for presence-observations
  #locate them to the closest cell coordinates
  ind_PO <- apply(species_data, 1, function(x) which.min((x[1]-env_data$coordinates[,1])^2 + 
                                                           (x[2]-env_data$coordinates[,2])^2))
  
  #define presence locations and covariates
  presence_coordinates <- env_data$coordinates[ind_PO,]
  presence_covariates <- env_data$covariates[ind_PO,]
  
  #define quadrature points:
  # fine scale point locations close to the observations
  # coarse scale point locations far from the observations
  
  #sample the closest ones to the presence observations
  #use radius for choosing the quadrature points
  radius <- 2
  env_ind_fine_2 <- apply(env_data$coordinates, 1, function(x) any(sqrt((x[1]-presence_coordinates[,1])^2 + 
    (x[2]-presence_coordinates[,2])^2)<=radius))
  
  quad_fine_coordinates <- env_data$coordinates[env_ind_fine_2==TRUE,]
  quad_fine_covariates <- env_data$covariates[env_ind_fine_2==TRUE,]
  
  #coarse scale point locations
  #load raster to get its resolution
  raster_cov_temp <- stack(paste(wd,"Data/Environment/Chelsa_SA.tif", sep = '/'))
  raster_res <- res(raster_cov_temp)/1000
  
  #create a grid for covariate data
  grid_sparse <- meshgrid(seq(raster_res[1]*scale_out/2, max(env_data$coordinates[,1]) - raster_res[1]*scale_out/2, raster_res[1]*scale_out),
                          seq(raster_res[2]*scale_out/2, max(env_data$coordinates[,2]) - raster_res[2]*scale_out/2, raster_res[2]*scale_out))
  
  #transform the grid into two vectors of coordinates
  env_coordinates_coarse <- cbind(array(grid_sparse$X), array(grid_sparse$Y))
  
  #transform the coordinates to original scale
  env_coordinates_coarse_t <- env_coordinates_coarse*1000
  env_coordinates_coarse_t <- cbind(env_coordinates_coarse_t[,1]+env_data$min_coordinates[1],
                                    env_coordinates_coarse_t[,2]+env_data$min_coordinates[2])
  
  #drop points which are outside of the species domain
  domain_temp <- stack(paste(wd,"/Data/Domains/", species, ".tif", sep = ''))
  ind_domain <- raster::extract(domain_temp[[2]],env_coordinates_coarse_t)
  env_coordinates_coarse_t <- env_coordinates_coarse_t[!is.na(ind_domain),]
  
  #assign covariate values to coarse grid cells
  coarse_scale_covariates <- raster::extract(raster_cov_temp,env_coordinates_coarse_t)
  
  #transform coordinates
  env_coordinates_coarse_t <- cbind(env_coordinates_coarse_t[,1]-env_data$min_coordinates[1],
                                    env_coordinates_coarse_t[,2]-env_data$min_coordinates[2])/1000

  #find rows where any covariate has NA value
  ind_NA <- apply(coarse_scale_covariates,1,anyNA)
  
  quad_coarse_coordinates <- env_coordinates_coarse_t[ind_NA==FALSE,]
  quad_coarse_covariates <- coarse_scale_covariates[ind_NA==FALSE,]
  
  #remove the cells which are too close to the fine scale quadrature points
  radius <- max(raster_res)*sqrt(2)*(scale_out/2+.5)
  
  ind_remove <- apply(quad_coarse_coordinates, 1,
                      function(x) any(sqrt((x[1] - quad_fine_coordinates[,1])^2 +
                                             (x[2] - quad_fine_coordinates[,2])^2)<radius))
  
  quad_coarse_coordinates <- quad_coarse_coordinates[ind_remove==FALSE,]
  quad_coarse_covariates <- quad_coarse_covariates[ind_remove==FALSE,]
  
  #combine all observations
  train_coordinates <- rbind(presence_coordinates,quad_fine_coordinates,quad_coarse_coordinates)
  train_covariates <- rbind(presence_covariates,quad_fine_covariates,quad_coarse_covariates)
  
  #standardize the covariates
  cov_mean <- apply(train_covariates,2,mean)
  cov_sd <- apply(train_covariates,2,sd)
  train_covariates_st <- apply(train_covariates, 2, function(x) (x-mean(x))/sd(x))
  
  #create a response variable
  response <-c(rep(1,nrow(presence_coordinates)),rep(0,nrow(quad_fine_coordinates)),
               rep(0,nrow(quad_coarse_coordinates)))
  
  #assign weights to the quadarature weights and the presence observations
  #weights in the presence locations is .5 since there is as well a quadrature point
  #in each presence location
  weights_presence <- rep(.5,nrow(presence_coordinates))
  
  #weights of the fine scale quadrature points
  weights_fine_quad <- rep(1, nrow(quad_fine_coordinates))
  
  #find quadrature points that are in the same locations
  ind_pres <- apply(presence_coordinates, 1, function(x) 
    which(x[1] == quad_fine_coordinates[,1] & x[2] == quad_fine_coordinates[,2]))
  
  #weights in those quadrature points are .5
  weights_fine_quad[ind_pres] <- .5
  
  #weights for sparse quadrature points are scale_out^2 or 1
  weight_coarse <- ifelse(weights_area, scale_out^2, 1)
  weights_coarse_quad <- rep(weight_coarse, nrow(quad_coarse_coordinates))

  #combine all weights into one vector
  weights <- c(weights_presence, weights_fine_quad, weights_coarse_quad)
  
  #define training data as a list
  training_data <- list(train_coordinates,train_covariates_st,cov_mean,cov_sd,env_data$min_coordinates,response,weights)
  names(training_data) <- c('coordinates', 'covariates', 'cov_mean', 'cov_sd', 'min_coordinates', 'response', 'weights')
  
  return(training_data)
}