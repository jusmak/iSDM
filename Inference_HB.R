Inference_HB <- function(wd, env_data, species_data) {
  
  
  set_scale = 100
  #create a grid for covariate data
  grid_sparse <- meshgrid(seq(0, max(env_data$coordinates[,1]), set_scale),
                          seq(0, max(env_data$coordinates[,2]), set_scale))
  
  env_coordinates_v2 <- cbind(array(grid_sparse$X), array(grid_sparse$Y))
  
  #get covariate values for the cells
  
  
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
  
  #extent of observations
  obs_extent <- c(min(presence_coordinates[,1]), max(presence_coordinates[,1]),
                  min(presence_coordinates[,2]), max(presence_coordinates[,2]))
  
  #fine scale point locations
  env_ind_fine <- 0
  for (i in 1:nrow(env_data$coordinates)){
    env_ind_fine[i] <- env_data$coordinates[i,1] >= obs_extent[1] & env_data$coordinates[i,1] <= obs_extent[2] &
    env_data$coordinates[i,2] >= obs_extent[3] & env_data$coordinates[i,2] <= obs_extent[4]
  }
  
  radius <- 5
  env_ind_fine <- 0
  for (i in 1:nrow(env_data$coordinates)){
    env_ind_fine[i] <- sum(sqrt((env_data$coordinates[i,1]-presence_coordinates[,1])^2 + 
      (env_data$coordinates[i,2]-presence_coordinates[,2])^2)<=radius)
  }
  
  #coarse scale point locations
  

  
  
  
}