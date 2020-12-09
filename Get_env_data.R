Get_env_data <- function(wd, species){
  
  #load environmental covariates as stacked rasters
  raster_cov_temp <- stack(paste(wd,"Data/Environment/Chelsa_SA.tif", sep = '/'))

  #load geographic extent
  domain_temp <- stack(paste(wd,"/Data/Domains/", species, ".tif", sep = ''))
  
  #crop covariates with the geographic extent
  raster_cov_temp <- crop(raster_cov_temp, domain_temp)
  
  #finalize cropping by setting all values outside of the extent to NA
  values(raster_cov_temp)[is.na(values(domain_temp)[,2]),] <- NA
  
  #create an empty matrix for storing values of the raster layers
  cov_matrix <- matrix(NA, nrow = length(values(raster_cov_temp[[1]])),
                       ncol = length(names(raster_cov_temp)))
  
  #store values
  for (i in 1:ncol(cov_matrix)) {
    cov_matrix[,i] <- values(raster_cov_temp[[i]])
  }
  
  #find rows where any covariate has NA and remove the row
  rows_na <- apply(cov_matrix,1,anyNA)
  cov_matrix <- cov_matrix[rows_na==FALSE,]
  
  #get coordinates of the cells
  xy_matrix <- coordinates(raster_cov_temp)[rows_na==FALSE,]
  
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
  covariates <- list(cov_matrix,xy_matrix,min_coord)
  names(covariates) <- c("covariates", "coordinates", "min_coordinates")
  
  return(covariates)
}