## iSDM: Compile data sets for multiple data types and conduct a joint inference with them

#---- Load packages ----#
suppressPackageStartupMessages({
  library(raster)
  library(pracma)
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
})



#---- Set working directory ----#
setwd('iSDM')
.wd = getwd()


#---- Initialize Environment ----#
.seed = 5326
set.seed(.seed)
options(mc.cores = parallel::detectCores())

#---- Set parameters for model runs ----#
#model_version
model_v = 1

#read in species list
species_list = read.csv('R/Species_list.csv', header = F)[[1]]

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


for (i in 1:length(species_list)) {
  .species = species_list[i]

  #get data
  source('R/Get_data.R')
  data = Get_data(.wd, .species, thinning, target_n_obs, weights_area, get_offsets)

  
  if (thinning_PA) {
    #thin species PA data
    #thinning and cleaning function
    source('R/occurrence.R')
    
    #save pa observations into a folder
    pa_temp = cbind(data$inventory_data$coordinates[,1]*1000 + data$training_data$min_coordinates[1],
                     data$inventory_data$coordinates[,2]*1000 + data$training_data$min_coordinates[2])
    pa_temp = cbind(pa_temp[,2],pa_temp[,1])
    colnames(pa_temp) = c('lat', 'lon')
    pa_folder = paste0('Data/PA_observations/Cleaned_observations/', .species, '.csv')
    write.csv(pa_temp, file = pa_folder)
    
    #define input data for thinning
    domain_temp = stack(paste0('Data/Domains/', .species, '.tif'))
    sp_pa_thin = cleanOcc(speciesCSV = pa_folder, env = domain_temp, doThin = thinning_PA)
    
    #derive points
    col_ind = match(c('pres.Longitude', 'pres.Latitude'),colnames(as.data.frame(sp_pa_thin)))
    sp_pa = as.data.frame(sp_pa_thin)[,col_ind]
    
    #scale the coordinates of the species occurrence records with the same coordinates
    sp_pa = cbind(sp_pa[,1]-data$training_data$min_coordinates[1],
                   sp_pa[,2]-data$training_data$min_coordinates[2])
    sp_pa = sp_pa/1000
    colnames(sp_pa) = c('X', 'Y')
    
    #find indexes of observations left after thinning
    ind_left = apply(sp_pa, 1, function(x) which.min((x[1]-data$inventory_data$coordinates[,1])^2 + 
                                                        (x[2]-data$inventory_data$coordinates[,2])^2))
    
    # subsample inventory data
    data$inventory_data$coordinates = data$inventory_data$coordinates[ind_left,]
    data$inventory_data$covariates = data$inventory_data$covariates[ind_left,]
    data$inventory_data$response = data$inventory_data$response[ind_left]
    for (j in 1:length(data$inventory_data$offset)) {
      data$inventory_data$offset[[j]] = data$inventory_data$offset[[j]][ind_left,]
    }
    
    data$pred_PA_data = data$inventory_data
  }

  #transfer presence-only data to correspond to the SPDE discratization
  source('R/Set_data_for_inference.R')
  data_inla = Set_data_spde(.wd, .species, thinning, target_n_obs, weights_area, data)
  
  #fit models
  source('R/RunInference_inla_non_spat.R')
  output_non_spat = RunInference(.wd, .species, data_inla, thinning, thinning_PA, weights_area, target_n_obs, model_v)
  
  #fit models
  source('R/RunInference_inla_spat.R')
  output_spat = RunInference(.wd, .species, data_inla, thinning, thinning_PA, weights_area, target_n_obs, model_v)
  
  #fit models
  source('R/RunInference_inla_spat_samp.R')
  output_spat_samp = RunInference(.wd, .species, data_inla, thinning, thinning_PA, weights_area, target_n_obs, model_v)
  
  #fit models for cv-samples
  source('R/RunInference_cv_inla_non_spat.R')
  output_cv_non_spat = RunInference_cv(.wd, .species, data_inla, thinning, thinning_PA, weights_area, target_n_obs, model_v)
  
  #fit models for cv-samples
  source(.wd, 'R/RunInference_cv_inla_spat.R')
  output_cv_spat = RunInference_cv(.wd, .species, data_inla, thinning, thinning_PA, weights_area, target_n_obs, model_v)
  
  #fit models for cv-samples
  source(.wd, 'R/RunInference_cv_inla_spat_samp.R')
  output_cv_spat_samp = RunInference_cv(.wd, .species, data_inla, thinning, thinning_PA, weights_area, target_n_obs, model_v)
  
  #fit models
  source('R/RunInference_inla_spat_rsr.R')
  output_spat_rsr = RunInference(.wd, .species, data_inla, thinning, thinning_PA, weights_area, target_n_obs, model_v)
  
  #fit models for cv-samples
  source('R/RunInference_cv_inla_spat_rsr.R')
  output_cv_spat_rsr = RunInference_cv(.wd, .species, data_inla, thinning, thinning_PA, weights_area, target_n_obs, model_v)
  
  #validate models of cv-samples
  source('R/Cross_validation.R')
  cross_validation = Validation_cv(.wd, .species, data, data_inla, thinning, thinning_PA, weights_area, target_n_obs, n_samples, model_v)
}

