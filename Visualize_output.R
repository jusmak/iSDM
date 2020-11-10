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

#Set SDM parameters
species = c("Loddigesia_mirabilis", "Ocreatus_underwoodii", "Oreotrochilus_leucopleurus")
thinning = FALSE
observations = "PO"
geographic_extent = "South-America"

#Loop over all species
for (i in 1:length(species)) {
  
  #Load data
  load(paste(.wd, '/', species[i], '_model_fits.RData', sep = ''))
  
  #Create predictions over the geographic extent
  #Get environmental prediction data
  source(paste(.wd, "R_code/Get_prediction_data.R", sep = '/'))
  pred_data <- Get_prediction_data(.wd, model_fits$training_data, species[i])
  
  #Compute predictive accuracy of the models
  #loo_cv,Tjur R for test and training data
  source(paste(.wd, "R_code/Predictive_tests.R", sep = '/'))
  pred_metric <- Predictive_test(.wd, model_fits, pred_data)
}


