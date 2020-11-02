## SDM ##
# Jussi Makinen 2020

# Define species, inference method and offset


doc <- 
'
Template

Usage:
script_template <out> [-t] [--seed=<seed>]
script_template (-h | --help)

Options:
-h --help     Show this screen.
-v --version     Show version.
-s --seed=<seed>  Random seed. Defaults to 5326 if not passed
-t --test         Indicates script is a test run, will not save output parameters or commit to git'

if(interactive()) {
  library(here)
  
  .wd <- 'Hummingbird_project'
  .script <- 'R_code/SDM.r' #Currently executing script
  .seed <- NULL
  .test <- TRUE
  rd <- here
  
  .outPF <- file.path(.wd,'Figures/myfig.png')
  
} else {
  library(docopt)
  library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  
  .wd <- getwd()
  .script <-  thisfile()
  .seed <- ag$seed
  .test <- as.logical(ag$test)
  rd <- is_rstudio_project$make_fix_file(.script)
  
  .outPF <- ag$out
}



#---- Initialize Environment ----#
.seed <- ifelse(is.null(.seed),5326,as.numeric(.seed))

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

# Set SDM parameters
species = "Oreotrochilus_leucopleurus"
offset = "none"
observations = "PO"
geographic_extent = "South-America"

source()

output <- RunInference(.wd, species, offset, observations, geographic_extent)
