Predictive_tests <- function(wd, model_fits, pred_PA_data) {

  #parallelize computation
  options(mc.cores = parallel::detectCores())
  
  #update here the real ind_fold
  #get occurrence points
  occ_points <- model_fits$training_data$coordinates[model_fits$training_data$response==1,]
  #transform the coordinates to the original scale
  occ_points <- cbind(occ_points[,1]*1000+model_fits$training_data$min_coordinates[1],
                      occ_points[,2]*1000+model_fits$training_data$min_coordinates[2])
  
  #get the projection of the env raster
  raster_cov_temp <- stack(paste(wd,"Data/Environment/Chelsa_SA.tif", sep = '/'))
  proj_raster <- projection(raster_cov_temp)
  
  #create a data frame for species names
  sp_name <- data.frame(rep(species, nrow(occ_points)))
  colnames(sp_name) <- 'sp'
  
  #transform occurrence points into a spatialpoints data frame
  occ_sp_points <- SpatialPointsDataFrame(coords = occ_points,
                                          data = sp_name,
                                          proj4string = crs(proj_raster))
  
  #do a clustering-based division to folds
  source(paste(wd, "R_code/spatialCV.r", sep = '/'))
  sp_stats <- list(species)
  names(sp_stats) <- 'species'
  n_fold <- spatialStratify(occ_sp_points, sp_stats, nfolds = -1, nsubclusters = -1)
  ind_fold_temp <- n_fold@data$folds
  #create a list of folds
  ind_fold <- list()
  for (i in 1:length(unique(ind_fold_temp))) {
    ind_fold[[i]] <- which(ind_fold_temp == i)
  }

  ##HB
  ##5-fold-cv
  source(paste(wd, "R_code/Inference_HB_cv_x_fold.R", sep = '/'))
  HB_fit_temp <- Inference_HB_cv_x_fold(wd, model_fits$training_data, model_fits$offset, ind_fold)

  return(HB_fit_temp)
  }


