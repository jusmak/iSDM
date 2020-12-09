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

#---- Set parameters for model runs ----#

#species
species_list = c("Loddigesia_mirabilis", "Ocreatus_underwoodii", "Oreotrochilus_leucopleurus")

#SDMs are fitted with each combination of offsets
#"none", "expert range map", "elevation" and "expert range map + elevation"

#option for spatially thin occurrence records
thinning = FALSE

#set the number of quadrature points in the coarse scale grid
targe_n_obs = 50000

#inference types
HB_inf = FALSE
glmnet_inf = FALSE
ppml_inf = TRUE

# RunInference_v1 takes in environmental and species data, runs the inference and saves the model output
source(paste(.wd, "R_code/RunInference_v1.R", sep = '/'))
for (i in 1:length(species_list)) {
  .species <- species_list[i]
  model_fits <- RunInference_v1(.wd, .species, thinning, targe_n_obs, HB_inf, glmnet_inf, ppml_inf)
  save(model_fits, file = paste(.wd, '/Model_fits/', .species, '_model_fits.RData', sep = ''))
}




