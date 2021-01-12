PP_training_data <- function(wd, env_data, species_data, scale_out, weights_area) {
  
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
  radius <- 5
  env_ind_fine_2 <- rep(NA,nrow(env_data$coordinates))
  for (i in 1:nrow(env_data$coordinates)){
    env_ind_fine_2[i] <- sum(sqrt((env_data$coordinates[i,1]-presence_coordinates[,1])^2 + 
      (env_data$coordinates[i,2]-presence_coordinates[,2])^2)<=radius)>0
  }
  quad_fine_coordinates <- env_data$coordinates[env_ind_fine_2==1,]
  quad_fine_covariates <- env_data$covariates[env_ind_fine_2==1,]
  
  
  #coarse scale point locations
  #load raster to get its resolution
  raster_cov_temp <- stack(paste(wd,"Data/Environment/Chelsa_SA.tif", sep = '/'))
  raster_res <- res(raster_cov_temp)/1000
  
  #create a grid for covariate data
  grid_sparse <- meshgrid(seq(raster_res[1]*scale_out/2, max(env_data$coordinates[,1]), raster_res[1]*scale_out),
                          seq(raster_res[2]*scale_out/2, max(env_data$coordinates[,2]), raster_res[2]*scale_out))
  
  #transform the grid into two vectors of coordinates
  env_coordinates_coarse <- cbind(array(grid_sparse$X), array(grid_sparse$Y))
  
  #transform the coordinates to original scale
  env_coordinates_coarse_t <- env_coordinates_coarse*1000
  env_coordinates_coarse_t <- cbind(env_coordinates_coarse_t[,1]+env_data$min_coordinates[1],
                                    env_coordinates_coarse_t[,2]+env_data$min_coordinates[2])
  
  #assign covariate values to coarse grid cells
  coarse_scale_covariates <- raster::extract(raster_cov_temp,env_coordinates_coarse_t)
  
  #find rows where any covariate has NA value
  ind_NA <- apply(coarse_scale_covariates,1,anyNA)
  
  quad_coarse_coordinates <- env_coordinates_coarse[ind_NA==FALSE,]
  quad_coarse_covariates <- coarse_scale_covariates[ind_NA==FALSE,]
  
  #remove the cells which are too close to the fine scale quadrature points
  radius <- max(raster_res)*(scale_out/2+.5)
  ind_remove <- rep(NA, nrow(quad_coarse_coordinates))
  for (i in 1:nrow(quad_coarse_coordinates)) {
    ind_remove[i] <- sum(sqrt((quad_coarse_coordinates[i,1] - quad_fine_coordinates[,1])^2 +
                                (quad_coarse_coordinates[i,2] - quad_fine_coordinates[,2])^2)<radius)>0
    
  }
  
  quad_coarse_coordinates <- quad_coarse_coordinates[ind_remove==0,]
  quad_coarse_covariates <- quad_coarse_covariates[ind_remove==0,]
  
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