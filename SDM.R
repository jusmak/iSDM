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

suppressPackageStartupMessages({
  library(raster)
  library(glmnet)
  library(pracma)
  library(ppmlasso)
  library(rstan)
  library(bayesplot)
  library(rstanarm)
  library(ggplot2)
  library(dismo)
})
#----  ----#

##SDMs are fitted with each combination of offsets 

# Set SDM parameters
species = c("Loddigesia_mirabilis", "Ocreatus_underwoodii", "Oreotrochilus_leucopleurus")
thinning = FALSE
observations = "PO"
geographic_extent = "South-America"
source(paste(.wd, "Rcode/RunInference_v1.R", sep = '/'))

for (i in length(species)) {
  model_fits <- RunInference_v1(.wd, species[i], observations, geographic_extent)
  save(model_fits, file = paste(.wd, '/', species[i], '_model_fits.RData', sep = ''))
}
t1 <- Sys.time()
