Get_env_data <- function(wd, geographic_extent){
  
  #load raster
  raster_cov_temp <- stack(paste(wd,"Data/Environment/Chelsa_Americas.tif", sep = '/'))
  
  #load geographic extent
  extent <- read.csv(paste(wd,"Data/Environment/Extent_geographic_area.csv", sep = '/'), sep = ';', header = TRUE)
  extent_temp <- as.double(unlist(extent[extent[,1]==geographic_extent,2:5]))
  
  #crop geographic area
  raster_cov_temp <- crop(raster_cov_temp, extent_temp)
  
  #name the raster layers
  names(raster_cov_temp) <- c("cov1", "cov2", "cov3", "cov4")
  
  #create an empty matrix for storing values of the raster layers
  cov_matrix <- matrix(NA, nrow = sum(!is.na(values(raster_cov_temp[[2]]))),
                       ncol = length(names(raster_cov_temp)))
  
  #store values
  for (i in 1:ncol(cov_matrix)) {
    cov_matrix[,i] <- values(raster_cov_temp[[i]])[!is.na(values(raster_cov_temp[[2]]))]
  }
  
  #find rows where any covariate has NA and remove the row
  rows_na <- apply(cov_matrix,1,anyNA)
  cov_matrix <- cov_matrix[rows_na==FALSE,]

  #get coordinates of the cells
  xy_matrix <- coordinates(raster_cov_temp)[!is.na(values(raster_cov_temp[[2]])),]
  xy_matrix <- xy_matrix[rows_na==FALSE,]

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