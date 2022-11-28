## Summaries for data and model outputs: species list with number of observations, average distance between observations,
## size of the species range from range maps, size of the predicted species range, species richness,
## predictive uncertainty

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


title_labels = c('AUC', 'TSS', 'Kappa', 'Specificity', 'Sensitivity', 'Log-lik', 'Tjur-R')
model_labels = c('PO', 'PA', 'Joint', 'PO', 'PA', 'PO+PA', 'PO+PA+samp.',
                 'PO + restr.', 'PA + restr.', 'PO+PA+restr.')
scenario_labels = c('Interpolation', 'Extrapolation')


#Make a list of species with the number of data points in both data sets
po_temp = matrix(NA, length(species_list), 2)
pa_temp = matrix(NA, length(species_list), 2)
source('R/Get_data.R')
source('R/Set_data_for_inference.R')

for (i in 1:length(species_list)) {
  .species = species_list[i]
  data = Get_data(.wd, .species, thinning, target_n_obs, weights_area)
  data_inla = Set_data_spde(.wd, .species, thinning, target_n_obs, weights_area, data)
  
  po_temp[i,1] = sum(data_inla$po_response==1)
  po_temp[i,2] = sum(data_inla$po_response==0)
  
  pa_temp[i,1] = sum(data_inla$pa_response==1)
  pa_temp[i,2] = sum(data_inla$pa_response==0)
}

species_n_data = list(PO = po_temp, PA = pa_temp)

save(species_n_data, file = 'Model_fits_v1/Summaries_cv_sp/Species_n_data.RData')


# average distance between data points in PO and PA data sets
dist_temp = matrix(NA, length(species_list), 1)

for (i in 1:length(species_list)) {
  .species = species_list[i]
  data = Get_data(.wd, .species, thinning, target_n_obs, weights_area)
  data_inla = Set_data_spde(.wd, .species, thinning, target_n_obs, weights_area, data)
  temp = 0
  
  for (j in 1:nrow(data_inla$pa_coordinates)) {
    for (k in 1:nrow(data_inla$po_coordinates[data_inla$po_response==1,])) {
      temp = temp + sqrt((data_inla$pa_coordinates[j,1] - data_inla$po_coordinates[k,1])^2 +
                            (data_inla$pa_coordinates[j,2] - data_inla$po_coordinates[k,2])^2)
    }
  }
  dist_temp[i] = temp/(nrow(data_inla$pa_coordinates)*nrow(data_inla$po_coordinates[data_inla$po_response==1,]))
}

save(dist_temp, file = 'Model_fits_v1/Summaries_cv_sp/Species_dist_data.RData')


# size of the species range (from expert range map)
range_size_temp = matrix(NA, length(species_list), 1)

for (i in 1:length(species_list)) {
  .species = species_list[i]
  
  #load geographic extent and expert range map
  domain_temp = stack(paste0(.wd,'/Data/Domains_new/', .species, '.tif'))
  expert_prior = raster(paste0(.wd, '/Data/Range_map_new/', .species, '.tif'))
  
  #crop the expert map
  expert_prior = crop(expert_prior, domain_temp)
  
  # get size of the range map
  range_size_temp[i] = sum(!is.na(values(expert_prior)))
}

save(range_size_temp, file = 'Model_fits_v1/Summaries_cv_sp/Species_range_data.RData')


#species richness from the range maps
# get extent of the study domain
extent_range = matrix(NA, length(species_list), 4)
for (i in 1:length(species_list)) {
  .species = species_list[i]
  #link range prediction to a large raster
  domain_temp = stack(paste0(.wd,'/Data/Domains_new/', .species, '.tif'))
  background = domain_temp
  extent_range[i,] = raster::extent(background)[1:4]
}
extent_max = c(min(extent_range[,1], na.rm = T), max(extent_range[,2], na.rm = T),
                min(extent_range[,3], na.rm = T), max(extent_range[,4], na.rm = T))

# distribution of the predicted species richness
species_rich_temp = raster(paste0(.wd, '/Data/Range_map_new/', species_list[1], '.tif'))
values(species_rich_temp) = NA
raster_temp = raster()
extent(raster_temp) = extent_max
species_rich = crop(species_rich_temp, raster_temp)

# get range map for each species
for (i in 1:length(species_list)) {
  .species = species_list[i]
  range_temp = raster(paste0(.wd, '/Data/Range_map_new/', species_list[i], '.tif'))
  range_temp = crop(range_temp, species_rich)
  species_rich = raster::mosaic(species_rich, range_temp, fun = 'sum')
}

save(species_rich, file = 'Model_fits_v1/Summaries_cv_sp/Species_richness_range_maps.RData')


# predictive uncertainty
model_titles = c('non_spat_po', 'non_spat_pa', 'non_spat_comb', 'spat_po', 'spat_pa', 'spat_comb', 'spat_comb_pref',
                  'spat_po_rsr', 'spat_pa_rsr', 'spat_comb_rsr')
model_ind = c(1:3,1:3,1,1:3)

sp_var = vector(mode = 'list', length = 8)
for (i in 1:length(sp_var)) {
  sp_var[[i]] = matrix(NA, length(species_list), length(model_titles))
}

for (i in 1:length(species_list)) {
  .species = species_list[i]
  print(2)
  print(i)
  for (j in 1:length(model_titles)) {
    
    if (j == 1) {
      load(file = paste0('Model_fits_v', model_v, '/INLA_non_spat_', .species, '_fits_CV_', thin_mark, thin_PA_mark, '_quad_n_',
                        target_n_obs, '_', weights_mark, '_inla.RData'))
    } else if (j == 4) {
      load(file = paste0('Model_fits_v', model_v, '/INLA_spat_', .species, '_fits_CV_', thin_mark, thin_PA_mark, '_quad_n_',
                        target_n_obs, '_', weights_mark, '_inla.RData'))
    } else if (j == 7) {
      load(file = paste0('Model_fits_v', model_v, '/INLA_spat_samp_', .species, '_fits_CV_', thin_mark, thin_PA_mark, '_quad_n_',
                        target_n_obs, '_', weights_mark, '_inla.RData'))
    } else if (j == 8) {
      load(file = paste0('Model_fits_v', model_v, '/INLA_spat_rsr_', .species, '_fits_CV_', thin_mark, thin_PA_mark, '_quad_n_',
                        target_n_obs, '_', weights_mark, '_inla.RData'))
    }      
    
    for (k in 1:length(sp_var)) {
      
      if (length(output_cv[[k]]$model_fits[[model_ind[j]]]) > 1) {
        #find the prediction cells
        ind_pred = inla.stack.index(output_cv[[k]]$stk_coll[[model_ind[j]]], 'pa_test')$data
        #collect their marginal standard deviation
        var_temp = output_cv[[k]]$model_fits[[model_ind[j]]]$summary.fitted.values[ind_pred,2]
        sp_var[[k]][i,j] = mean(var_temp)
      }
    }
  }
}

names(sp_var) = c('interpolation1', 'interpolation2', 'interpolation3', 'interpolation4',
                   'extrapolation1', 'extrapolation2', 'extrapolation3', 'extrapolation4')
save(sp_var, file = 'Model_fits_v1/Summaries_cv_sp/Species_var_pred.RData')



## Get model parameter estimates
model_titles = c('non_spat_po', 'non_spat_pa', 'non_spat_comb', 'spat_po', 'spat_pa', 'spat_comb', 'spat_comb_pref',
                  'spat_po_rsr', 'spat_pa_rsr', 'spat_comb_rsr')
model_ind = c(1:3,1:3,1,1:3)

range_est = matrix(NA, length(species_list)*n_samples, length(model_titles))
sd_est = matrix(NA, length(species_list)*n_samples, length(model_titles))
range_samp_est = matrix(NA, length(species_list)*n_samples, 1)
sd_samp_est = matrix(NA, length(species_list)*n_samples, 1)

for (i in 1:length(species_list)) {
  .species = species_list[i]
  for (j in 4:length(model_titles)) {
    
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
    
    #find the prediction cells
    if (length(output$model_fits[[model_ind[j]]]) > 1) {
      range_est[((i-1)*n_samples+1):(i*n_samples),j] = inla.rmarginal(n_samples, output$model_fits[[model_ind[j]]]$marginals.hyperpar[[1]])
      sd_est[((i-1)*n_samples+1):(i*n_samples),j] = inla.rmarginal(n_samples, output$model_fits[[model_ind[j]]]$marginals.hyperpar[[2]])
      
      if (j == 7) {
        range_samp_est[((i-1)*n_samples+1):(i*n_samples)] = inla.rmarginal(n_samples, output$model_fits[[model_ind[j]]]$marginals.hyperpar[[3]])
        sd_samp_est[((i-1)*n_samples+1):(i*n_samples)] = inla.rmarginal(n_samples, output$model_fits[[model_ind[j]]]$marginals.hyperpar[[4]])
      }
    }
  }
}

estimates = list(range = range_est, sd = sd_est, range_samp = range_samp_est, sd_samp = sd_samp_est)
save(estimates, file = 'Model_fits_v1/Summaries_cv_sp/Estimates.RData')


# size of the predicted range
pred_range_size_temp = matrix(NA, length(species_list), 3)
source('R/Get_data.R')
source('R/Set_data_for_inference.R')

model_titles = c('non_spat_po', 'non_spat_pa', 'non_spat_comb', 'spat_po', 'spat_pa', 'spat_comb', 'spat_comb_pref',
                  'spat_po_rsr', 'spat_pa_rsr', 'spat_comb_rsr')
model_ind = c(1:3,1:3,1,1:3)

return_thres = function(response, pred_temp) {
  my_roc = suppressMessages(pROC::roc(response = response, predictor = pred_temp))
  return(my_roc$thresholds[which.max(my_roc$sensitivities+my_roc$specificities)[1]])
}

# get extent of the study domain
extent_range = matrix(NA, length(species_list), 4)
for (i in 1:length(species_list)) {
  .species = species_list[i]
  #link range prediction to a large raster
  domain_temp = stack(paste0('Data/Domains_new/', .species, '.tif'))
  background = domain_temp
  extent_range[i,] = raster::extent(background)[1:4]
}
extent_max = c(min(extent_range[,1], na.rm = T), max(extent_range[,2], na.rm = T),
                min(extent_range[,3], na.rm = T), max(extent_range[,4], na.rm = T))

# distribution of the predicted species richness
species_rich_temp = raster(paste0('Data/Range_map_new/', species_list[1], '.tif'))
values(species_rich_temp) = NA
raster_temp = raster()
extent(raster_temp) = extent_max
species_rich = crop(species_rich_temp, raster_temp)
species_rich_list = list(species_rich, species_rich, species_rich, species_rich,
                          species_rich, species_rich, species_rich, species_rich,
                          species_rich, species_rich)
pred_range_size_temp = matrix(NA, nrow = length(species_list), ncol = length(model_titles))
species_rich_list_2 = list(species_rich, species_rich, species_rich, species_rich,
                            species_rich, species_rich, species_rich, species_rich,
                            species_rich, species_rich)
pred_range_size_temp_2 = matrix(NA, nrow = length(species_list), ncol = length(model_titles))

for (i in 1:length(species_list)) {
  .species = species_list[i]
  print(1)
  print(i)
  data = Get_data(.wd, .species, thinning, target_n_obs, weights_area)
  data_inla = Set_data_spde(.wd, .species, thinning, target_n_obs, weights_area, data)
  
  #predictive map
  #set prediction cells into a matrix
  pred_matrix = data$env_data$covariates
  
  #take the second order effect
  pred_matrix = cbind(pred_matrix, pred_matrix^2)
  
  #scale the covariate values according to the training locations
  for (j in 1:ncol(pred_matrix)) {
    pred_matrix[,j] = (pred_matrix[,j]-data_inla$cov_mean[j])/data_inla$cov_sd[j]
  }
  
  #add vector of ones for a constant
  pred_matrix = cbind(rep(1,nrow(pred_matrix)),pred_matrix)
  
  #get presence-absence and presence-only cells for setting threshold
  pred_matrix_pa = cbind(rep(1,length(data_inla$pa_response)), data_inla$pa_covariates)
  
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
      inla_temp = output$model_fits[[model_ind[j]]]
      
      #predict across the study area
      if (j == 1 | j == 2 | j == 4 | j == 5 | j == 8 | j == 9) {
        beta_mean = inla_temp$summary.fixed$mean
      } else {
        beta_mean = inla_temp$summary.fixed$mean[-2]
      }
      
      #covariate effect prediction
      pred_mean = pred_matrix%*%as.matrix(beta_mean)
      
      #spatial random effect prediction
      if (j == 5 | j == 9) {
        proj = inla.mesh.projector(data_inla$spde_mesh, data$env_data$coordinates)
        sre_mean = inla.mesh.project(proj, output$model_fits[[model_ind[j]]]$summary.random$pa_spde$mean)
        pred_mean = pred_mean + sre_mean
      } else if (j == 4 | j == 6 | j == 7 | j == 8 | j == 10) {
        proj = inla.mesh.projector(data_inla$spde_mesh, data$env_data$coordinates)
        sre_mean = inla.mesh.project(proj, output$model_fits[[model_ind[j]]]$summary.random$po_spde$mean)
        pred_mean = pred_mean + sre_mean
      }
      
      #get a threshold for the predictions by using PA training points
      #covariate effect prediction
      pred_mean_thres = pred_matrix_pa%*%as.matrix(beta_mean)
      
      #spatial random effect prediction
      if (j == 5 | j == 9) {
        proj = inla.mesh.projector(data_inla$spde_mesh, data_inla$pa_coordinates)
        sre_mean = inla.mesh.project(proj, output$model_fits[[model_ind[j]]]$summary.random$pa_spde$mean)
        pred_mean_thres = pred_mean_thres + sre_mean
      } else if (j == 4 | j == 6 | j == 7 | j == 8 | j == 10) {
        proj = inla.mesh.projector(data_inla$spde_mesh, data_inla$pa_coordinates)
        sre_mean = inla.mesh.project(proj, output$model_fits[[model_ind[j]]]$summary.random$po_spde$mean)
        pred_mean_thres = pred_mean_thres + sre_mean
      }
      
      thres = return_thres(data_inla$pa_response, as.vector(pred_mean_thres))
      
      pred_cat = ifelse(pred_mean >= thres, 1, 0)
      
      #range prediction
      pred_range_size_temp[i,j] = sum(pred_cat, na.rm = T)
      
      #link range prediction to a large raster
      domain_temp = stack(paste0('Data/Domains_new/', .species, '.tif'))
      background = domain_temp
      values(background)[-data$env_data$ind_na] = NA
      background[data$env_data$ind_na] = pred_cat
      
      species_rich_list[[j]] = raster::mosaic(species_rich_list[[j]], background, fun = 'sum')
      
      if (j == 2 | j == 5 | j == 6 | j == 7 | j == 9 | j == 10) {
        #range prediction
        pred_range_size_temp_2[i,j] = sum(1-exp(-exp(pred_mean)) >= .5, na.rm = T)
        
        #link range prediction to a large raster
        domain_temp = stack(paste0('Data/Domains_new/', .species, '.tif'))
        background = domain_temp
        values(background)[-data$env_data$ind_na] = NA
        background[data$env_data$ind_na] = 1-exp(-exp(pred_mean))
        
        species_rich_list_2[[j]] = raster::mosaic(species_rich_list_2[[j]], background, fun = 'sum')
      }
    }
  }
}

pred_range = list(range = pred_range_size_temp, richness_map = species_rich_list, range_2 = pred_range_size_temp_2, richness_map_2 = species_rich_list_2)

save(pred_range, file = 'Model_fits_v1/Summaries_cv_sp/Species_range_pred.RData')


# proportion of predicted species areal densities falling into the range
source('R/Get_data.R')
source('R/Set_data_for_inference.R')

model_titles = c('non_spat_po', 'non_spat_pa', 'non_spat_comb', 'spat_po', 'spat_pa', 'spat_comb', 'spat_comb_pref',
                  'spat_po_rsr', 'spat_pa_rsr', 'spat_comb_rsr')
model_ind = c(1:3,1:3,1,1:3)

dens_range_map = matrix(NA, length(species_list), length(model_titles))


for (i in 1:length(species_list)) {
  .species = species_list[i]
  print(1)
  print(i)
  data = Get_data(.wd, .species, thinning, target_n_obs, weights_area)
  data_inla = Set_data_spde(.wd, .species, thinning, target_n_obs, weights_area, data)
  
  #predictive map
  #set prediction cells into a matrix
  pred_matrix = data$env_data$covariates
  
  #take the second order effect
  pred_matrix = cbind(pred_matrix, pred_matrix^2)
  
  #scale the covariate values according to the training locations
  for (j in 1:ncol(pred_matrix)) {
    pred_matrix[,j] = (pred_matrix[,j]-data_inla$cov_mean[j])/data_inla$cov_sd[j]
  }
  
  #add vector of ones for a constant
  pred_matrix = cbind(rep(1,nrow(pred_matrix)),pred_matrix)
  
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
      inla_temp = output$model_fits[[model_ind[j]]]
      
      #predict across the study area
      if (j == 1 | j == 2 | j == 4 | j == 5 | j == 8 | j == 9) {
        beta_mean = inla_temp$summary.fixed$mean
      } else {
        beta_mean = inla_temp$summary.fixed$mean[-2]
      }
      
      #covariate effect prediction
      pred_mean = pred_matrix%*%as.matrix(beta_mean)
      
      #spatial random effect prediction
      if (j == 5 | j == 9) {
        proj = inla.mesh.projector(data_inla$spde_mesh, data$env_data$coordinates)
        sre_mean = inla.mesh.project(proj, output$model_fits[[model_ind[j]]]$summary.random$pa_spde$mean)
        pred_mean = pred_mean + sre_mean
      } else if (j == 4 | j == 6 | j == 7 | j == 8 | j == 10) {
        proj = inla.mesh.projector(data_inla$spde_mesh, data$env_data$coordinates)
        sre_mean = inla.mesh.project(proj, output$model_fits[[model_ind[j]]]$summary.random$po_spde$mean)
        pred_mean = pred_mean + sre_mean
      }
      
      #get species range map
      range_temp = raster(paste0('Data/Range_map_new/', .species, '.tif'))
      #link range map to domain
      domain_temp = stack(paste0('Data/Domains_new/', .species, '.tif'))
      background = domain_temp
      values(background)[-data$env_data$ind_na] = NA
      background[data$env_data$ind_na] = exp(pred_mean)
      range_temp = crop(range_temp, background)
      range_temp[-data$env_data$ind_na] = NA
      dens_range_map[i,j] = sum(values(background)[!is.na(values(range_temp))])/sum(values(background), na.rm = T)
      
    }
  }
}

dens_range = list(range = dens_range_map)

save(dens_range, file = 'Model_fits_v1/Summaries_cv_sp/Dens_range_map.RData')


# locations of all PO observations
PO_temp = matrix(NA, nrow = 1, ncol = 2)

for (i in 1:length(species_list)) {
  .species = species_list[i]
  data = Get_data(.wd, .species, thinning, target_n_obs, weights_area)
  data_inla = Set_data_spde(.wd, .species, thinning, target_n_obs, weights_area, data)
  coord_temp = cbind(data_inla$po_coordinates[data_inla$po_response==1,1]*1000 + data$env_data$min_coordinates[1],
                     data_inla$po_coordinates[data_inla$po_response==1,2]*1000 + data$env_data$min_coordinates[2])
  PO_temp = rbind(PO_temp, coord_temp)
}
PO_temp = PO_temp[-1,]

save(PO_temp, file = 'Model_fits_v1/Summaries_cv_sp/Species_po_data.RData')






