## Compute the predictive accuracy of range predictions by using range maps as the validation data
## Two prediction scales: 1km and 50km


#---- Load packages ----#
suppressPackageStartupMessages({
  library(raster)
  library(ggplot2)
  library(caret)
  library(pROC)
  library(flexclust)
  library(irr)
  library(spdep)
  library(INLA)
  library(gridExtra)
  library(gtable)
  library(RColorBrewer)
  library(wesanderson)
  library(tidyverse)
  library(rgdal)
  library(R.utils)
  library(maptools)
})



#---- Set working directory ----#
setwd('iSDM')
.wd = getwd()


#---- Initialize Environment ----#
.seed = 5326
set.seed(.seed)
options(mc.cores = parallel::detectCores())

#---- Set parameters for model runs ----#
#model version
model_v = 1

#read in species list
species_list = read.csv('R/Species_list_cross_validation.csv', header = F)[[1]]

#option for spatially thin PO records
thinning = FALSE
thinning_PA = FALSE

#set the number of quadrature points in the coarse scale grid
target_n_obs = 20000

#size of predictive marginal distributions
n_samples = 500

#get offset values
get_offsets = FALSE

#whether weights are assigned according to the area related to a quadrature point
weights_area = TRUE

#create indices
thin_mark = ifelse(thinning, 'thin', 'not_thin')
thin_PA_mark =  ifelse(thinning_PA, '_thin_PA', '_not_thin_PA')
weights_mark = ifelse(weights_area, 'wa', 'not_wa')


source('R/Get_data.R')
source('R/Set_data_for_inference.R')

model_labels = c('PO', 'PA', 'Joint', 'PO', 'PA', 'PO+PA', 'PO+PA+samp.',
                 'PO + restr.', 'PA + restr.', 'PO+PA+restr.')

scenario_labels = c('Interpolation', 'Extrapolation')
model_ind <- c(1:3,1:3,1,1:3)

return_thres <- function(response, pred_temp) {
  my_roc <- suppressMessages(pROC::roc(response = response, predictor = pred_temp))
  return(my_roc$thresholds[which.max(my_roc$sensitivities+my_roc$specificities)[1]])
}


# 50km
range_accuracy <- vector(mode = 'list', length(species_list))

metric_labels = c('AUC', 'Kappa', 'TSS', 'Specificity', 'Sensitivity')

for (i in 1:length(species_list)) {
  accuracy_table <- matrix(NA, length(metric_labels), length(model_titles))
  rownames(accuracy_table) = metric_labels
  colnames(accuracy_table) = model_titles
  
  .species <- species_list[i]
  print(i)
  data <- Get_data(.wd, .species, thinning, target_n_obs, weights_area)
  data_inla <- Set_data_spde(.wd, .species, thinning, target_n_obs, weights_area, data)
  
  #predictive map
  #set prediction cells into a matrix
  pred_matrix <- data$env_data$covariates
  
  #take the second order effect
  pred_matrix <- cbind(pred_matrix, pred_matrix^2)
  
  #scale the covariate values according to the training locations
  for (j in 1:ncol(pred_matrix)) {
    pred_matrix[,j] <- (pred_matrix[,j]-data_inla$cov_mean[j])/data_inla$cov_sd[j]
  }
  
  #add vector of ones for a constant
  pred_matrix <- cbind(rep(1,nrow(pred_matrix)),pred_matrix)
  
  #get presence-absence and presence-only cells for setting threshold
  pred_matrix_pa <- cbind(rep(1,length(data_inla$pa_response)), data_inla$pa_covariates)
  
  #get species range map
  range_temp <- raster(paste0('Data/Range_map_new/', .species, '.tif'))
  raster::values(range_temp)[!is.na(raster::values(range_temp))] = 1
  raster::values(range_temp)[is.na(raster::values(range_temp))] = 0
  raster::values(range_temp)[-data$env_data$ind_na] <- NA
  
  #upscale range map
  range_temp_agg = aggregate(range_temp, fact=50, fun=max)
  range_map_all <- raster::values(range_temp_agg)
  
  #keep only 0s and 1s
  range_map <- range_map_all[-which(is.na(range_map_all))]
  
  
  
  #loop over different models
  for (j in 1:length(model_titles)) {
    
    if (j == 1) {
      load(file = paste0('Model_fits_v', model_v, '/INLA_non_spat_', .species, '_fits_', thin_mark, thin_PA_mark, '_quad_n_',
                         target_n_obs, '_', weights_mark, '_inla.RData'))
    } else if (j == 4) {
      load(file = paste0('Model_fits_v', model_v, '/INLA_spat_', .species, '_fits_', thin_mark, thin_PA_mark, '_quad_n_',
                         target_n_obs, '_', weights_mark, '_inla.RData'))
    } else if (j == 7) {
      load(file = paste0('Model_fits_v', model_v, '/INLA_spat_samp_', .species, '_fits_', thin_mark, thin_PA_mark, '_quad_n_',
                         target_n_obs, '_', weights_mark, '_inla.RData'))
    } else if (j == 8) {
      load(file = paste0('Model_fits_v', model_v, '/INLA_spat_rsr_', .species, '_fits_', thin_mark, thin_PA_mark, '_quad_n_',
                         target_n_obs, '_', weights_mark, '_inla.RData'))
    }
    
    #pick a model
    if (length(output$model_fits[[model_ind[j]]]) > 1) {
      inla_temp <- output$model_fits[[model_ind[j]]]
      
      #predict across the study area
      if (j == 1 | j == 2 | j == 4 | j == 5 | j == 8 | j == 9) {
        beta_mean <- inla_temp$summary.fixed$mean
      } else {
        beta_mean <- inla_temp$summary.fixed$mean[-2]
      }
      
      #covariate effect prediction
      pred_mean <- pred_matrix%*%as.matrix(beta_mean)
      
      #spatial random effect prediction
      if (j == 5 | j == 9) {
        proj <- inla.mesh.projector(data_inla$spde_mesh, data$env_data$coordinates)
        sre_mean <- inla.mesh.project(proj, output$model_fits[[model_ind[j]]]$summary.random$pa_spde$mean)
        pred_mean <- pred_mean + sre_mean
      } else if (j == 4 | j == 6 | j == 7 | j == 8 | j == 10) {
        proj <- inla.mesh.projector(data_inla$spde_mesh, data$env_data$coordinates)
        sre_mean <- inla.mesh.project(proj, output$model_fits[[model_ind[j]]]$summary.random$po_spde$mean)
        pred_mean <- pred_mean + sre_mean
      }
      
      #use range map as a base for upscaling predictions
      raster::values(range_temp)[data$env_data$ind_na] <- pred_mean
      pred_mean_agg = raster::aggregate(range_temp, fact=50, fun=mean)
      pred_mean_agg = raster::values(pred_mean_agg)
      pred_mean_agg <- pred_mean_agg[-which(is.na(range_map_all))]
      
      #continuous predictor
      my_roc = suppressMessages(pROC::roc(response = range_map, predictor = pred_mean_agg))
      accuracy_table[1,j] = my_roc$auc
      
      #get a threshold for the predictions by using PA training points
      #covariate effect prediction
      pred_mean_thres <- pred_matrix_pa%*%as.matrix(beta_mean)
      
      #spatial random effect prediction
      if (j == 5 | j == 9) {
        proj <- inla.mesh.projector(data_inla$spde_mesh, data_inla$pa_coordinates)
        sre_mean <- inla.mesh.project(proj, output$model_fits[[model_ind[j]]]$summary.random$pa_spde$mean)
        pred_mean_thres <- pred_mean_thres + sre_mean
      } else if (j == 4 | j == 6 | j == 7 | j == 8 | j == 10) {
        proj <- inla.mesh.projector(data_inla$spde_mesh, data_inla$pa_coordinates)
        sre_mean <- inla.mesh.project(proj, output$model_fits[[model_ind[j]]]$summary.random$po_spde$mean)
        pred_mean_thres <- pred_mean_thres + sre_mean
      }
      
      thres <- return_thres(data_inla$pa_response, as.vector(pred_mean_thres))
      pred_cat <- ifelse(pred_mean >= thres, 1, 0)
      
      raster::values(range_temp)[data$env_data$ind_na] <- pred_cat
      pred_cat_agg = raster::aggregate(range_temp, fact=50, fun=max)
      pred_cat_agg = raster::values(pred_cat_agg)
      pred_cat_agg = pred_cat_agg[-which(is.na(range_map_all))]
      
      
      #categorical predictor
      accuracy_table[2,j] <- kappa2(cbind(range_map,pred_cat_agg))$value
      ind_thres <- which.min(abs(my_roc$thresholds-thres)) 
      accuracy_table[3,j] <- my_roc$sensitivities[ind_thres]+my_roc$specificities[ind_thres]-1
      accuracy_table[4,j] <- my_roc$specificities[ind_thres]
      accuracy_table[5,j] <- my_roc$sensitivities[ind_thres]
      
    }
  }
  range_accuracy[[i]] = accuracy_table
}

save(range_accuracy, file = 'Model_fits_v1/Summaries_cv_sp/Range_map_accuracy_50km.RData')


# 1km
range_accuracy <- vector(mode = 'list', length(species_list))

metric_labels = c('AUC', 'Kappa', 'TSS', 'Specificity', 'Sensitivity')

for (i in 1:length(species_list)) {
  accuracy_table <- matrix(NA, length(metric_labels), length(model_titles))
  rownames(accuracy_table) = metric_labels
  colnames(accuracy_table) = model_titles
  
  .species <- species_list[i]
  print(i)
  data <- Get_data(.wd, .species, thinning, target_n_obs, weights_area)
  data_inla <- Set_data_spde(.wd, .species, thinning, target_n_obs, weights_area, data)
  
  #predictive map
  #set prediction cells into a matrix
  pred_matrix <- data$env_data$covariates
  
  #take the second order effect
  pred_matrix <- cbind(pred_matrix, pred_matrix^2)
  
  #scale the covariate values according to the training locations
  for (j in 1:ncol(pred_matrix)) {
    pred_matrix[,j] <- (pred_matrix[,j]-data_inla$cov_mean[j])/data_inla$cov_sd[j]
  }
  
  #add vector of ones for a constant
  pred_matrix <- cbind(rep(1,nrow(pred_matrix)),pred_matrix)
  
  #get presence-absence and presence-only cells for setting threshold
  pred_matrix_pa <- cbind(rep(1,length(data_inla$pa_response)), data_inla$pa_covariates)
  
  #get species range map
  range_temp <- raster(paste(.wd, '/Data/Range_map_new/', .species, '.tif', sep = ''))
  raster::values(range_temp)[is.na(raster::values(range_temp))] = 0
  raster::values(range_temp)[-data$env_data$ind_na] <- NA
  range_map <- raster::values(range_temp)
  #keep only 0s and 1s
  range_map <- range_map[-which(is.na(range_map))]
  
  #loop over different models
  for (j in 1:length(model_titles)) {
    
    if (j == 1) {
      load(file = paste0('Model_fits_v', model_v, '/INLA_non_spat_', .species, '_fits_', thin_mark, thin_PA_mark, '_quad_n_',
                         target_n_obs, '_', weights_mark, '_inla.RData'))
    } else if (j == 4) {
      load(file = paste0('Model_fits_v', model_v, '/INLA_spat_', .species, '_fits_', thin_mark, thin_PA_mark, '_quad_n_',
                         target_n_obs, '_', weights_mark, '_inla.RData'))
    } else if (j == 7) {
      load(file = paste0('Model_fits_v', model_v, '/INLA_spat_samp_', .species, '_fits_', thin_mark, thin_PA_mark, '_quad_n_',
                         target_n_obs, '_', weights_mark, '_inla.RData'))
    } else if (j == 8) {
      load(file = paste0('Model_fits_v', model_v, '/INLA_spat_rsr_', .species, '_fits_', thin_mark, thin_PA_mark, '_quad_n_',
                         target_n_obs, '_', weights_mark, '_inla.RData'))
    }
    
    #pick a model
    if (length(output$model_fits[[model_ind[j]]]) > 1) {
      inla_temp <- output$model_fits[[model_ind[j]]]
      
      #predict across the study area
      if (j == 1 | j == 2 | j == 4 | j == 5 | j == 8 | j == 9) {
        beta_mean <- inla_temp$summary.fixed$mean
      } else {
        beta_mean <- inla_temp$summary.fixed$mean[-2]
      }
      
      #covariate effect prediction
      pred_mean <- pred_matrix%*%as.matrix(beta_mean)
      
      #spatial random effect prediction
      if (j == 5 | j == 9) {
        proj <- inla.mesh.projector(data_inla$spde_mesh, data$env_data$coordinates)
        sre_mean <- inla.mesh.project(proj, output$model_fits[[model_ind[j]]]$summary.random$pa_spde$mean)
        pred_mean <- pred_mean + sre_mean
      } else if (j == 4 | j == 6 | j == 7 | j == 8 | j == 10) {
        proj <- inla.mesh.projector(data_inla$spde_mesh, data$env_data$coordinates)
        sre_mean <- inla.mesh.project(proj, output$model_fits[[model_ind[j]]]$summary.random$po_spde$mean)
        pred_mean <- pred_mean + sre_mean
      }
      
      #continuous predictor
      my_roc = suppressMessages(pROC::roc(response = range_map, predictor = pred_mean))
      accuracy_table[1,j] = my_roc$auc
      
      #get a threshold for the predictions by using PA training points
      #covariate effect prediction
      pred_mean_thres <- pred_matrix_pa%*%as.matrix(beta_mean)
      
      #spatial random effect prediction
      if (j == 5 | j == 9) {
        proj <- inla.mesh.projector(data_inla$spde_mesh, data_inla$pa_coordinates)
        sre_mean <- inla.mesh.project(proj, output$model_fits[[model_ind[j]]]$summary.random$pa_spde$mean)
        pred_mean_thres <- pred_mean_thres + sre_mean
      } else if (j == 4 | j == 6 | j == 7 | j == 8 | j == 10) {
        proj <- inla.mesh.projector(data_inla$spde_mesh, data_inla$pa_coordinates)
        sre_mean <- inla.mesh.project(proj, output$model_fits[[model_ind[j]]]$summary.random$po_spde$mean)
        pred_mean_thres <- pred_mean_thres + sre_mean
      }
      
      thres <- return_thres(data_inla$pa_response, as.vector(pred_mean_thres))
      pred_cat <- ifelse(pred_mean >= thres, 1, 0)
      
      #categorical predictor
      accuracy_table[2,j] <- kappa2(cbind(range_map,pred_cat))$value
      ind_thres <- which.min(abs(my_roc$thresholds-thres)) 
      accuracy_table[3,j] <- my_roc$sensitivities[ind_thres]+my_roc$specificities[ind_thres]-1
      accuracy_table[4,j] <- my_roc$specificities[ind_thres]
      accuracy_table[5,j] <- my_roc$sensitivities[ind_thres]
      
    }
  }
  range_accuracy[[i]] = accuracy_table
}

save(range_accuracy, file = 'Model_fits_v1/Summaries_cv_sp/Range_map_accuracy_1km.RData')


