## Visualize output ##
# Jussi Makinen 2020

#Define species, inference method and offset

if(interactive()) {
  .wd <- 'C:/Users/OMISTAJA/Documents/Hummingbird_project'
} else {
  .wd <- '/gpfs/loomis/pi/jetz/jm3669/Hummingbird_project'
}

.script <- 'SDM.r' #Currently executing script

#---- Initialize Environment ----#
.seed <- 5326
set.seed(.seed)
if(interactive()) {
  lib_path <- .libPaths()[1]
} else {
  lib_path <- "/gpfs/loomis/pi/jetz/jm3669/R/3.6/"
}

suppressPackageStartupMessages({
  library(raster, lib.loc=lib_path)
  library(glmnet, lib.loc=lib_path)
  library(pracma, lib.loc=lib_path)
  library(ppmlasso, lib.loc=lib_path)
  library(rstan, lib.loc=lib_path)
  library(bayesplot, lib.loc=lib_path)
  library(rstanarm, lib.loc=lib_path)
  library(ggplot2, lib.loc=lib_path)
  library(dismo, lib.loc=lib_path)
  library(dplyr, lib.loc=lib_path)
  library(loo, lib.loc=lib_path)
  library(caret, lib.loc=lib_path)
  library(pROC, lib.loc=lib_path)
  library(flexclust, lib.loc=lib_path)
  library(RColorBrewer, lib.loc=lib_path)
})

#Set SDM parameters
species_list = c("Loddigesia_mirabilis", "Ocreatus_underwoodii", "Oreotrochilus_leucopleurus") #
observations = "PO"
geographic_extent = "South-America"
thinning = FALSE
n_samp = 500
targe_n_obs = 50000

#Loop over all species
for (i in 1:length(species_list)) {
  
  .species <- species_list[i]
  
  #set scale out according to the target_n_obs
  domain_temp <- stack(paste(.wd,"/Data/Domains/", .species, ".tif", sep = ''))
  scale_out <- round(sqrt(sum(!is.na(values(domain_temp[[1]])))/targe_n_obs))
  
  #Load data
  load(paste(.wd, '/', .species, '_model_fits.RData', sep = ''))
  print(.species)
  
  #Create predictions over the geographic extent
  #Get environmental prediction data
  source(paste(.wd, "R_code/Get_prediction_data.R", sep = '/'))
  pred_data <- Get_prediction_data(.wd, model_fits$training_data, .species, scale_out)
  print("Pred_data derived")
  
  #Get validation data
  source(paste(.wd, "R_code/Get_validation_data.R", sep = '/'))
  pred_PA_data <- Get_validation_data(.wd, model_fits$training_data, model_fits$offset, .species)
  print("Pred_PA_data derived")
  
  #Compute predictive accuracy of the models
  #cv for training data and predictive density for validation data
  #source(paste(.wd, "R_code/Predictive_tests.R", sep = '/'))
  #pred_metric <- Predictive_tests(.wd, model_fits, .species, pred_PA_data)
  
  #save(pred_metric, file = paste(.wd, '/', .species, '_validation.RData', sep = ''))
  
  #Visualizations
  #source(paste(.wd, "R_code/Prediction_visualization_HB.R", sep = '/'))
  #Prediction_visualization(.wd, model_fits, pred_data, n_samp, .species, scale_out)

  source(paste(.wd, "R_code/Prediction_visualization_glmnet.R", sep = '/'))
  Prediction_visualization(.wd, model_fits, pred_data, n_samp, .species, scale_out)
  
}


