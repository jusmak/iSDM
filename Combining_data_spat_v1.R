## Combining data test #1 v1 ##
# Jussi Makinen 2021

#---- Load packages ----#
if(interactive()) {
  lib_path <- .libPaths()[1]
} else {
  lib_path <- .libPaths()[1]
#  lib_path <- "/gpfs/loomis/pi/jetz/jm3669/R/3.6/"
}

suppressPackageStartupMessages({
  library(raster, lib.loc=lib_path)
  library(glmnet, lib.loc=lib_path)
  library(pracma, lib.loc=lib_path)
  library(rstan, lib.loc=lib_path)
  library(rstanarm, lib.loc=lib_path)
  library(ggplot2, lib.loc=lib_path)
  library(dismo, lib.loc=lib_path)
  library(dplyr, lib.loc=lib_path)
  library(loo, lib.loc=lib_path)
  library(caret, lib.loc=lib_path)
  library(pROC, lib.loc=lib_path)
  library(flexclust, lib.loc=lib_path)
  library(irr, lib.loc=lib_path)
})


#---- Set working directory ----#
if(interactive()) {
  .wd <- 'C:/Users/OMISTAJA/OneDrive - Yale University/Hummingbird_project'
} else {
  .wd <- '/gpfs/loomis/pi/jetz/jm3669/Hummingbird_project'
}

#---- Initialize Environment ----#
.seed <- 5326
set.seed(.seed)
options(mc.cores = parallel::detectCores())

#---- Set parameters for model runs ----#

#species
species_list = c('Adelomyia_melanogenys',
                 'Aglaeactis_cupripennis',
                 'Amazilia_amazilia',
                 'Amazilia_franciae',
                 'Amazilia_saucerrottei',
                 'Amazilia_tzacatl',
                 'Anthracothorax_nigricollis',
                 'Boissonneaua_flavescens',
                 'Chaetocercus_mulsant',
                 'Chalybura_buffonii',
                 'Chlorostilbon_mellisugus',
                 'Coeligena_coeligena',
                 'Coeligena_torquata',
                 'Colibri_coruscans',
                 'Colibri_thalassinus',
                 'Doryfera_ludovicae',
                 'Ensifera_ensifera',
                 'Eutoxeres_aquila',
                 'Florisuga_mellivora',
                 'Glaucis_hirsutus',
                 'Haplophaedia_aureliae',
                 'Heliangelus_exortis',
                 'Heliomaster_longirostris',
                 'Lafresnaya_lafresnayi',
                 'Metallura_tyrianthina',
                 'Ocreatus_underwoodii',
                 'Phaethornis_guy',
                 'Phaethornis_striigularis',
                 'Phaethornis_syrmatophorus',
                 'Ramphomicron_microrhynchum',
                 'Thalurania_furcata')

#SDMs are fitted without offsets

#option for spatially thin PO records
thinning = FALSE
thinning_PA_opt = FALSE

#set the number of quadrature points in the coarse scale grid
target_n_obs = 10000

#in hierarchical Bayesian inference:
# full data: 5000 samples with 2000 burn-in
# cross-validation folds: 1200 samples with 400 burn-in
n_samples_1 = 5000
burn_in_1 = 2000
n_samples_2 = 1500
burn_in_2 = 500

#whether weights are assigned according to the area related to a quadrature point
weights_area_opt = TRUE


for (t in 1:length(thinning_PA_opt)) {
  thinning_PA = thinning_PA_opt[t]
  
  for (w in 1:length(weights_area_opt)) {
    weights_area = weights_area_opt[w]
    
    for (i in 1:length(species_list)) {
      .species <- species_list[i]
    
      #get data
      source(paste(.wd, "R_code/Get_data.R", sep = '/'))
      data <- Get_data(.wd, .species, thinning, target_n_obs, weights_area)

      if (thinning_PA) {
        #thin species PA data
        #thinning and cleaning function
        source(paste(.wd, "R_code/occurrence.R", sep = '/'))
        
        #save pa observations into a folder
        pa_temp <- cbind(data$inventory_data$coordinates[,1]*1000 + data$training_data$min_coordinates[1],
                         data$inventory_data$coordinates[,2]*1000 + data$training_data$min_coordinates[2])
        pa_temp <- cbind(pa_temp[,2],pa_temp[,1])
        colnames(pa_temp) <- c('lat', 'lon')
        pa_folder <- paste(.wd, "/Data/PA_observations/Cleaned_observations/", .species, ".csv", sep = '')
        write.csv(pa_temp, file = pa_folder)
        
        #define input data for thinning
        domain_temp <- stack(paste(.wd,"/Data/Domains/", .species, ".tif", sep = ''))
        sp_pa_thin <- cleanOcc(speciesCSV = pa_folder, env = domain_temp, doThin = thinning_PA)
        
        #derive points
        col_ind <- match(c("pres.Longitude", "pres.Latitude"),colnames(as.data.frame(sp_pa_thin)))
        sp_pa <- as.data.frame(sp_pa_thin)[,col_ind]
        
        #scale the coordinates of the species occurrence records with the same coordinates
        sp_pa <- cbind(sp_pa[,1]-data$training_data$min_coordinates[1],
                       sp_pa[,2]-data$training_data$min_coordinates[2])
        sp_pa <- sp_pa/1000
        colnames(sp_pa) <- c('X', 'Y')
        
        #find indexes of observations left after thinning
        ind_left <- apply(sp_pa, 1, function(x) which.min((x[1]-data$inventory_data$coordinates[,1])^2 + 
                                                            (x[2]-data$inventory_data$coordinates[,2])^2))
        
        # subsample inventory data
        data$inventory_data$coordinates <- data$inventory_data$coordinates[ind_left,]
        data$inventory_data$covariates <- data$inventory_data$covariates[ind_left,]
        data$inventory_data$response <- data$inventory_data$response[ind_left]
        for (j in 1:length(data$inventory_data$offset)) {
          data$inventory_data$offset[[j]] <- data$inventory_data$offset[[j]][ind_left,]
        }
        
        data$pred_PA_data <- data$inventory_data
      }
      
      #set inventory observations into training and testing sets
      #get occurrence points
      occ_points <- data$inventory_data$coordinates
      #transform the coordinates to the original scale
      occ_points <- cbind(occ_points[,1]*1000+data$training_data$min_coordinates[1],
                          occ_points[,2]*1000+data$training_data$min_coordinates[2])
      
      #create a data frame for species names
      sp_name <- data.frame(rep(.species, nrow(occ_points)))
      colnames(sp_name) <- 'sp'
      
      #transform occurrence points into a spatialpoints data frame
      occ_sp_points <- SpatialPointsDataFrame(coords = occ_points,
                                              data = sp_name,
                                              proj4string = crs(data$proj_raster))
      
      #do a clustering-based division to folds
      source(paste(.wd, "R_code/spatialCV.r", sep = '/'))
      sp_stats <- list(.species)
      names(sp_stats) <- 'species'
      n_fold <- spatialStratify(occ_sp_points, sp_stats, nfolds = 4, nsubclusters = 3)
      ind_fold_temp <- n_fold@data$folds
      #create a list of folds
      ind_fold <- list()
      for (j in 1:length(unique(ind_fold_temp))) {
        ind_fold[[j]] <- which(ind_fold_temp == j)
      }
      
      # add flag if a subsample has less than 3 presence records
      pres_temp <- 0
      for (j in 1:length(ind_fold)) {
        pres_temp[j] <- sum(data$inventory_data$response[ind_fold[[j]]])
      }
      presence_flag <- ifelse(any(pres_temp < 3), 1, 0)
      
      #fit models
      source(paste(.wd, "R_code/RunInference_v2.R", sep = '/'))
      model_fits <- RunInference_v2(.wd, .species, data, thinning, thinning_PA, weights_area, target_n_obs, n_samples_1, burn_in_1)
      
      #compute predictive accuracy of the models
      #source(paste(.wd, "R_code/Predictive_tests_test_3.R", sep = '/'))
      #pred_metric <- Predictive_tests(.wd, model_fits, .species, data, thinning, thinning_PA,
      #                                weights_area, target_n_obs, presence_flag, n_samples_2, burn_in_2)
    }
  }
}

