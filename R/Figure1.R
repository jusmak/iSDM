## Manuscript: Figure 1

#Function for finding the threshold for plot
#log-intensity value where sensitivity drops below sens (manually set sensitivity)
return_thres = function(response, pred_temp) {
  my_roc = suppressMessages(pROC::roc(response = response, predictor = pred_temp))
  return(my_roc$thresholds[which.max(my_roc$sensitivities+my_roc$specificities)[1]])
}

#---- Load packages ----#
suppressPackageStartupMessages({
  library(raster)
  library(ggplot2)
  library(pracma)
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


#set the color scale
col.brk=data.frame(cols=c(grey(seq(.8,.2,length=24)),colorRampPalette(c('steelblue4', 'steelblue1','gold','gold3', 'red1', 'red4'))(50)))
#make sure the top color is red
col.brk$cols[(nrow(col.brk)-1)]='#8B0000'

env_temp = raster('Study_area/Study_area_laea.tif')
#aggregate to 25x25km resolution
data('wrld_simpl')
wrld_simpl = spTransform(wrld_simpl, crs(env_temp))

#load the range map based validation
load(file = 'Model_fits_v1/Summaries_cv_sp/Range_map_accuracy_50km.RData')
load(file = 'Model_fits_v1/Summaries_cv_sp/Dens_range_map.RData')


#Choose two species
#Coeligena lutetiae (7)
#Ensifera ensifera (10)

#figures for the species (using PO, PA and PO+PA+samp models)
# 1. range and observations
# 2. range predictions
# 3. interpolation results

for (i in c(7,10)) {
  .species = species_list[i]
  
  ## Load data and model fits
  source('R/Get_data.R')
  data = Get_data(.wd, .species, thinning, target_n_obs, weights_area, get_offsets)
  
  # Load data set specific to INLA run
  source('R/Set_data_for_inference.R')
  data_inla = Set_data_spde(.wd, .species, thinning, target_n_obs, weights_area, data)
  
  #### Figure ####
  obs_points_1 = data.frame(cbind(data_inla$po_coordinates[data_inla$po_response==1,1] * 1000 + data$env_data$min_coordinates[1],
                                   data_inla$po_coordinates[data_inla$po_response==1,2] * 1000 + data$env_data$min_coordinates[2]))
  obs_points_2 = data.frame(cbind(data_inla$pa_coordinates[,1] * 1000 + data$env_data$min_coordinates[1],
                                   data_inla$pa_coordinates[,2] * 1000 + data$env_data$min_coordinates[2]))

  #load expert range map
  expert_temp = readOGR(paste0('Data/Range_shp/', .species, '.shp'))
  env_temp[is.na(env_temp[])] = 0 
  
  #transform the projection of the range map
  expert_temp = spTransform(expert_temp, crs(env_temp))
  
  #Fig1: range and po observations
  #plot the po+range map
  if (i == 10) {
    x1 = c(-86, -55)
    y1 = c(-20, 15)
  } else if (i == 7) {
    x1 = c(-86, -55)
    y1 = c(-8, 15)
  }

  d = data.frame(lon=x1, lat=y1)
  coordinates(d) = c('lon', 'lat')
  proj4string(d) = CRS('+init=epsg:4326') # WGS 84
  CRS.new = crs(env_temp)
  d = spTransform(d, CRS.new)
  coord_temp = d@coords
  
  #crop the expert_temp and env_temp
  raster_temp = raster()
  extent(raster_temp) = c(coord_temp[c(1,2),1], coord_temp[c(1,2),2])
  env_temp_2 = crop(env_temp, raster_temp)
  
  if (i == 10) {
    pdf(paste0('Figures/Fig1/Fig1_sp', num2str(i, fmt = 0), '.pdf'), width=5, height=5)
  } else if (i == 7) {
    pdf(paste0('Figures/Fig1/Fig1_sp', num2str(i, fmt = 0), '.pdf'), width=5, height=4)
  }
  
  par(pty = 's', mai = c(.6,.6,.1,.1))
  
  #plot expert range map
  plot(env_temp_2, legend = F, axes = F, box = F, col = c('light blue', brewer.pal(n = 8, name = 'Set2')[7]))
  plot(wrld_simpl, add = T)
    
  #add training points
  points(obs_points_1, pch = 20, cex = 1)
  
  plot(expert_temp, col = '#11440050', add = T)

  if (i == 10) {
    x2 = c(-85, -60)
    y2 = c(-20, 10)
  } else if (i == 7) {
    x2 = c(-85, -60)
    y2 = c(-5, 15)
  }

  d = data.frame(lon=x2, lat=y2)
  coordinates(d) = c('lon', 'lat')
  proj4string(d) = CRS('+init=epsg:4326') # WGS 84
  CRS.new = crs(env_temp)
  d = spTransform(d, CRS.new)
  coord_temp_axis = d@coords
  
  axis(1, at=coord_temp_axis[,1], labels=c('85\u00B0 W', '60\u00B0 W'), pos = coord_temp[1,2]-1e5, cex.axis = 1.2,  tck = -.01)
  axis(2, at=coord_temp_axis[,2], labels=c(paste0(abs(y2[1]),'\u00B0 S'), paste0(y2[2], '\u00B0 N')), pos = coord_temp[1,1]-1e5, cex.axis = 1.2, tck = -.01)
  
  #add number of data points in the figure
  text(x=1.3e6, y=0e6, labels=paste0('n(PO)=', sum(data_inla$po_response)), cex=1.2)
  text(x=1.3e6, y=-.25e6, labels=paste0('prev.(PA)=', num2str(sum(data_inla$pa_response)/length(data_inla$pa_response),fmt=2)), cex=1.2)
  
  dev.off()
  

  #Fig2: range and pa observations
  if (i == 10) {
    x1 = c(-86, -55)
    y1 = c(-20, 15)
  } else if (i == 7) {
    x1 = c(-86, -55)
    y1 = c(-8, 15)
  }
  
  d = data.frame(lon=x1, lat=y1)
  coordinates(d) = c('lon', 'lat')
  proj4string(d) = CRS('+init=epsg:4326') # WGS 84
  CRS.new = crs(env_temp)
  d = spTransform(d, CRS.new)
  coord_temp = d@coords
  
  #crop the expert_temp and env_temp
  raster_temp = raster()
  extent(raster_temp) = c(coord_temp[c(1,2),1], coord_temp[c(1,2),2])
  env_temp_2 = crop(env_temp, raster_temp)
  
  pdf(paste0('Figures/Fig1/Fig2_sp', num2str(i, fmt = 0), '.pdf'), width=5, height=4)
  
  par(pty = 's', mai = c(.6,.6,.1,.1))
  
  #plot expert range map
  plot(env_temp_2, legend = F, axes = F, box = F, col = c('light blue', brewer.pal(n = 8, name = 'Set2')[7]))
  plot(wrld_simpl, add = T)
  
  #add training points
  points(obs_points_2[data_inla$pa_response==0,], pch = 16, cex = .5, col =  brewer.pal(n = 8, name = 'Dark2')[8])
  points(obs_points_2[data_inla$pa_response==1,], pch = 16, cex = .5)
  
  plot(expert_temp, col = '#11440050', add = T)
  
  dev.off()
  
  
  
  #Fig3: Log-intensity prediction
  if (i == 10) {
    x1 = c(-86, -55)
    y1 = c(-26, 15)
  } else if (i == 7) {
    x1 = c(-86, -55)
    y1 = c(-15, 15)
  }
  
  d = data.frame(lon=x1, lat=y1)
  coordinates(d) = c('lon', 'lat')
  proj4string(d) = CRS('+init=epsg:4326') # WGS 84
  CRS.new = crs(env_temp)
  d = spTransform(d, CRS.new)
  coord_temp = d@coords
  
  #crop the expert_temp and env_temp
  raster_temp = raster()
  extent(raster_temp) = c(coord_temp[c(1,2),1], coord_temp[c(1,2),2])
  env_temp_2 = crop(env_temp, raster_temp)
  
  pred_matrix = data$env_data$covariates
  
  #take the second order effect
  pred_matrix = cbind(pred_matrix, pred_matrix^2)
  #scale the covariate values according to the training locations
  for (j in 1:ncol(pred_matrix)) {
    pred_matrix[,j] = (pred_matrix[,j]-data_inla$cov_mean[j])/data_inla$cov_sd[j]
  }
  
  pred_matrix = cbind(rep(1,nrow(pred_matrix)),pred_matrix)
  pred_matrix_thres = cbind(rep(1,length(data_inla$pa_response)), data_inla$pa_covariates)
  
  #estimate/load model
  source('R/RunInference_inla_spat.R')
  output_spat = RunInference(.wd, .species, data_inla, thinning, thinning_PA, weights_area, target_n_obs, model_v)
  
  source('R/RunInference_inla_spat_samp.R')
  output_spat_samp = RunInference(.wd, .species, data_inla, thinning, thinning_PA, weights_area, target_n_obs, model_v)
  
  background = stack(paste0('Data/Domains_new/', .species, '.tif'))
  background_temp = aggregate(background, fact = 25)
  
  for (k in c(4,5,7)) {
    
    if (k==4 | k == 5) {
        output = output_spat
      } else if (k == 7) {
        output = output_spat_samp
    }
    
    model_ind = c(1:3,1:3,1,1:3)[k]
    
    #load rasters for visualizing the predictions
    background = background_temp
    values(background)[-data$env_data$ind_na] = NA
    
    #pick a model
    inla_temp = output$model_fits[[model_ind]]
    
    #predict across the study area
    if (k == 1 | k == 2 | k == 4 | k == 5 | k == 8 | k == 9) {
      beta_mean = inla_temp$summary.fixed$mean
    } else {
      beta_mean = inla_temp$summary.fixed$mean[-2]
    }
    
    #covariate effect prediction
    pred_mean = pred_matrix%*%as.matrix(beta_mean)
    
    #spatial random effect prediction
    if (k > 3) {
      proj = inla.mesh.projector(data_inla$spde_mesh, data$env_data$coordinates)
      if (k == 5 | k == 9) {
        sre_mean = inla.mesh.project(proj, output$model_fits[[model_ind]]$summary.random$pa_spde$mean)
      } else {
        sre_mean = inla.mesh.project(proj, output$model_fits[[model_ind]]$summary.random$po_spde$mean)
      }      
      pred_mean = pred_mean + sre_mean
    }
    
    pred_mean[pred_mean<quantile(pred_mean, .1, na.rm = T)] = quantile(pred_mean,.1, na.rm = T)
    
    #presence-absence map with observations (test set)
    values(background)[data$env_data$ind_na] = pred_mean
    
    if (i == 10) {
      pdf(paste0('Figures/Fig1/Fig3_', num2str(i, fmt = 0), '_model_', num2str(k, fmt = 0), '.pdf'), width=5, height=6)
    } else if (i == 7){
      pdf(paste0('Figures/Fig1/Fig3_', num2str(i, fmt = 0), '_model_', num2str(k, fmt = 0), '.pdf'), width=5, height=4)
    }
    
    par(pty = 's', mai = c(.2,.1,.1,.1))
    plot(env_temp_2, legend = F, axes = F, box = F, col = c('light blue', brewer.pal(n = 8, name = 'Set2')[7]))
    plot(background, legend = F, axes = F, box = F, add = T, col = col.brk$cols)
    plot(wrld_simpl, add = T)
    plot(background, legend.only=TRUE, horizontal = F, col = col.brk$cols, legend.width=.5, legend.shrink=0.5,
         smallplot=c(.18,.2, .15, .25), axis.args = list(cex.axis=.8))
    
    
    rect(.5e6,-.3e6,1.9e6,.7e6, col= rgb(1,1,1,alpha=0.35), border = NA)
    
    #add AUC
    load(file = paste0('Model_fits_v', model_v, '/INLA_', .species, '_validation_', thin_mark, thin_PA_mark, '_quad_n_',
                      target_n_obs, '_', weights_mark, '_inla.RData'))
    
    #extrapolation
    scenario_ext = grep('Ext', names(output_validation_cv))
    text(x=1.1e6, y=.5e6, labels=bquote(paste('AUC'['CV']*'=', .(round(mean(unlist(lapply(lapply(output_validation_cv, function(x)  unlist(x[k])), function(x) x[[1]]))[scenario_ext]),2)))), cex=1.2)
    
    #add AUC of range prediction
    #load the proportions
    if (i == 10){
      text(x=1e6, y=.2e6, labels=bquote(paste('AUC'['range']*'=', .(round(range_accuracy[[i]][1,k],2)))), cex=1.2)
    } else if (i == 7) {
      text(x=1.1e6, y=.2e6, labels=bquote(paste('AUC'['range']*'=', .(round(range_accuracy[[i]][1,k],2)))), cex=1.2)
    }
    
    
    #add proportion of species population inside the range map
    text(x=1.1e6, y=-.1e6, labels=paste0('P.=', round(dens_range$range[i,k],2)), cex=1.2)
    
    
    dev.off()
  }
}


