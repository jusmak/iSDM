## SDM ##
# Jussi Makinen 2020

#---- Load packages ----#
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
})

#---- Set working directory ----#
if(interactive()) {
  .wd <- 'C:/Users/OMISTAJA/Documents/Hummingbird_project'
} else {
  .wd <- '/gpfs/loomis/pi/jetz/jm3669/Hummingbird_project'
}

#---- Initialize Environment ----#
.seed <- 5326
set.seed(.seed)
options(mc.cores = parallel::detectCores())

#---- Set parameters for model runs ----#

#species
species_list = c("Loddigesia_mirabilis", "Ocreatus_underwoodii", "Oreotrochilus_leucopleurus")

#SDMs are fitted with each combination of offsets
#"none", "expert range map", "elevation" and "expert range map + elevation"

#option for spatially thin occurrence records
thinning = FALSE

#set the number of quadrature points in the coarse scale grid
target_n_obs = 50000

#set the number of samples taken randomly from the posterior distribution
#for computing the posterior predictive distribution
n_samp = 500

#in hierarchical Bayesian inference:
# full data: 5000 samples with 2000 burn-in
# cross-validation folds: 1200 samples with 400 burn-in
# there use posterior mean values 

#the SDM workflow is created for two model comparisons:
#1. Hierarchical Bayesian vs. glmnet inference
#2. Presence-only vs. combined likelihoods, both fitted with ppmlasso
m_category = 2

#model fits for each species represent either category 1 or 2

if (m_category == 1) {
  #inference types
  HB_inf = TRUE
  glmnet_inf = TRUE
  ppml_inf = FALSE
  #combine likelihoods for presence-only and presence-absence data
  lik_comb = FALSE
} else {
  #inference types
  HB_inf = FALSE
  glmnet_inf = FALSE
  ppml_inf = TRUE
  #combine likelihoods for presence-only and presence-absence data
  lik_comb = TRUE
}


# RunInference_v1 takes in environmental and species data, runs the inference and saves the model output
for (i in 1:length(species_list)) {
  .species <- species_list[i]
  
  #get data
  source(paste(.wd, "R_code/Get_data.R", sep = '/'))
  data <- Get_data(.wd, .species, thinning, target_n_obs)
  
  #fit models
  source(paste(.wd, "R_code/RunInference_v1.R", sep = '/'))
  model_fits <- RunInference_v1(.wd, .species, data, HB_inf, glmnet_inf, ppml_inf, lik, thinning, target_n_obs, m_category)
  
  #compute predictive accuracy of the models
  source(paste(.wd, "R_code/Predictive_tests.R", sep = '/'))
  pred_metric <- Predictive_tests(.wd, model_fits, .species, data$training_data, data$offset_full,
                                  data$pred_PA_data, data$proj_raster, HB_inf, glmnet_inf, ppml_inf,
                                  thinning, target_n_obs, m_category, n_samp)
  
  #visualizations
  if (HB_inf) {
    source(paste(.wd, "R_code/Prediction_visualization_HB.R", sep = '/'))
    Prediction_visualization(.wd, model_fits, data$training_data, data$offset_full, data$pred_data, n_samp, .species, data$scale_out, pred_metric)
  }
  
  if (glmnet_inf) {
    source(paste(.wd, "R_code/Prediction_visualization_glmnet.R", sep = '/'))
    Prediction_visualization(.wd, model_fits, data$training_data, data$offset_full, data$pred_data, n_samp, .species, data$scale_out, pred_metric)
  }
#  
#  if (ppml_inf) {
#    source(paste(.wd, "R_code/Prediction_visualization_ppml.R", sep = '/'))
#    Prediction_visualization(.wd, model_fits, data$pred_data, n_samp, .species, data$scale_out)
#  }  
  
  
}




