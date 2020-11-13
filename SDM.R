## SDM ##
# Jussi Makinen 2020

# Define species, inference method and offset

if(interactive()) {
  .wd <- 'C:/Users/OMISTAJA/Documents/Hummingbird_project'
} else {
  .wd <- '/gpfs/loomis/pi/jetz/jm3669/Hummingbird_project'
}

.script <- 'SDM.r' #Currently executing script

#---- Initialize Environment ----#
.seed <- 5326
set.seed(.seed)
t0 <- Sys.time()
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

#----  ----#

##SDMs are fitted with each combination of offsets 

# Set SDM parameters
species = c("Loddigesia_mirabilis", "Ocreatus_underwoodii", "Oreotrochilus_leucopleurus")
observations = "PO"
geographic_extent = "South-America"
thinning = FALSE
scale_out = 20   #set the factor by which the spatial grid is coarsened
source(paste(.wd, "R_code/RunInference_v1.R", sep = '/'))

for (i in length(species)) {
  model_fits <- RunInference_v1(.wd, species[i], observations, geographic_extent, thinning, scale_out)
  save(model_fits, file = paste(.wd, '/', species[i], '_model_fits.RData', sep = ''))
}
t1 <- Sys.time()


