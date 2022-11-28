## Manuscript: Figure 2

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

model_ind_all <- c(4:7,10)

#get predicted range map size
load(file = 'Model_fits_v1/Summaries_all_sp/Species_range_pred_50km.RData')
load(file = 'Model_fits_v1/Summaries_all_sp/Species_richness_range_maps_50km.RData')
load(file = 'Model_fits_v1/Summaries_all_sp/Species_range_data.RData')



##Species richness prediction

#change the ocean cells to NA
env_raster <- stack('Study_area/Study_area_laea.tif')
#crop raster with the richness map
env_raster_c <- raster::crop(env_raster[[1]], pred_range$richness_map[[1]])
#find cells which have NA value
ind_na <- which(is.na(raster::values(env_raster_c)))
raster::values(env_raster_c) = ifelse(!is.na(raster::values(env_raster_c)), 1, 0)

#set the color scale for absolute values
col.brk=data.frame(cols=c(grey(seq(.8,.2,length=24)),colorRampPalette(c('steelblue4', 'steelblue1','gold','gold3', 'red1', 'red4'))(50)))
#make sure the top color is red
col.brk$cols[(nrow(col.brk)-1)]='#8B0000'

#set the color scale for relative values
col.brk.rel=data.frame(cols=c(colorRampPalette(c('steelblue4', 'steelblue1'))(15), grey(.9), colorRampPalette(c('gold3', 'red1'))(15)))


## SET 1 ## Absolute values
#coordinate limits  
xl_min <- -1.3e6
xl_max <- 1e6
yl_min <- -.5e6
yl_max <- 2e6

r = raster(paste0(.wd, '/Data/Study_area/Study_area_laea.tif'))
data("wrld_simpl")
wrld_simpl <- spTransform(wrld_simpl, raster::crs(r))

x <- c(-85,-68)
y <- c(-7.5,10)
d <- data.frame(lon=x, lat=y)
coordinates(d) <- c("lon", "lat")
proj4string(d) <- CRS("+init=epsg:4326") # WGS 84
CRS.new <- raster::crs(r)
coord_trans <- spTransform(d, CRS.new)
coord_temp <- coord_trans@coords

r_agg = crop(r, pred_range$richness_map[[1]])
r_agg = raster::aggregate(r_agg, fact = 50)
ind_na_agg = which(is.na(raster::values(r_agg)))
env_raster_c = raster::aggregate(env_raster_c, fact = 50)

#with axes - scale to max
for (i in 1:length(model_ind_all)) {
  temp_1 <- pred_range$richness_map[[model_ind_all[i]]]
  temp_1[ind_na_agg] = NA
  max_z = max(raster::values(temp_1), na.rm = T)
  pdf(paste0('Figures/Fig3/Abs_model_', num2str(model_ind_all[i],0), '_axes.pdf'), width = 5, height = 5)
  par(pty = 's', mai = c(.6,.2,.2,.2))
  plot(env_raster_c, axes = F, box = F, legend = F, col = c('light blue', brewer.pal(n = 8, name = "Set2")[7]),
       ylim = c(yl_min,yl_max), xlim = c(xl_min, xl_max))
  plot(temp_1, axes = F, box = F, legend = F, col = col.brk$cols, zlim = c(0, max_z),  ylim = c(yl_min,yl_max), xlim = c(xl_min, xl_max),
       legend.width=.35, legend.shrink=0.5, axis.args=list(cex.axis=.8), add = T)
  axis(1, at=coord_temp[,1], labels=c('85\u00B0 W', '68\u00B0 W'), pos = yl_min, cex.axis = .8, tck = -.01)
  axis(2, at=coord_temp[,2], labels=c('7.5\u00B0 S', '10\u00B0 N'), pos = xl_min, cex.axis = .8, tck = -.01)
  dev.off()
}

#without axes - scale to max
for (i in 1:length(model_ind_all)) {
  temp_1 <- pred_range$richness_map[[model_ind_all[i]]]
  temp_1[ind_na_agg] = NA
  max_z = max(raster::values(temp_1), na.rm = T)
  pdf(paste0('Figures/Fig3/Abs_model_', num2str(model_ind_all[i],0), '_no_axes.pdf'), width = 5, height = 5)
  par(pty = 's', mai = c(.6,.2,.2,.2))
  plot(env_raster_c, axes = F, box = F, legend = F, col = c('light blue', brewer.pal(n = 8, name = "Set2")[7]),
       ylim = c(yl_min,yl_max), xlim = c(xl_min, xl_max))
  if (i==4) {
    plot(temp_1, axes = F, box = F, col = col.brk$cols, zlim = c(0, max_z),  ylim = c(yl_min,yl_max), xlim = c(xl_min, xl_max),
         legend.width=.35, legend.shrink=0.5, axis.args=list(cex.axis=.8, at = c(0,max_z), labels = c('min', 'max')), add = T)
  } else {
    plot(temp_1, axes = F, box = F, legend = F, col = col.brk$cols, zlim = c(0, max_z),  ylim = c(yl_min,yl_max), xlim = c(xl_min, xl_max),
         legend.width=.35, legend.shrink=0.5, axis.args=list(cex.axis=.8), add = T)
  }
  dev.off()
}


#plot the species richness from range maps
xl_min <- -1.3e6
xl_max <- 1e6
yl_min <- -.5e6
yl_max <- 2e6

temp_1 <- species_rich
temp_1[ind_na_agg] = NA
max_z = max(raster::values(temp_1), na.rm = T)
pdf('Figures/Fig3/Sp_richness_range_map.pdf', width = 5, height = 6)
par(pty = 's', mai = c(.2,.6,.2,.2))
plot(env_raster_c, axes = F, box = F, legend = F, col = c('light blue', brewer.pal(n = 8, name = "Set2")[7]),
     ylim = c(yl_min,yl_max), xlim = c(xl_min, xl_max))
plot(temp_1, axes = F, box = F, legend = F, col = col.brk$cols, zlim = c(0, max_z),  ylim = c(yl_min,yl_max), xlim = c(xl_min, xl_max),
     legend.width=.35, legend.shrink=0.5, axis.args=list(cex.axis=.8), add = T)
axis(1, at=coord_temp[,1], labels=c('85\u00B0 W', '68\u00B0 W'), pos = yl_min, cex.axis = .8, tck = -.01)
axis(2, at=coord_temp[,2], labels=c('7.5\u00B0 S', '10\u00B0 N'), pos = xl_min, cex.axis = .8, tck = -.01)
dev.off()


#plot the species observation density map
load(file = 'Model_fits_v1/Summaries_all_sp/Species_po_data.RData')

temp <- pred_range$richness_map[[model_ind_all[1]]]
temp[ind_na_agg] = NA
temp_obs = temp
counts = log(table(raster::cellFromXY(temp_obs,PO_temp)))
values(temp_obs)[as.numeric(names(counts))] = counts
values(temp_obs)[-as.numeric(names(counts))] = 0
temp_obs[ind_na_agg] = NA

xl_min <- -1.3e6
xl_max <- 1e6
yl_min <- -.5e6
yl_max <- 2e6


pdf('Figures/Fig3/Obs_density.pdf', width = 5, height = 5)
par(pty = 's', mai = c(.2,.6,.2,.2))
plot(env_raster_c, axes = F, box = F, legend = F, col = c('light blue', brewer.pal(n = 8, name = "Set2")[7]),
     ylim = c(yl_min,yl_max), xlim = c(xl_min, xl_max))
plot(temp_obs, axes = F, box = F, col = brewer.pal(n = 9, name = "YlGnBu")[3:9], ylim = c(yl_min,yl_max), xlim = c(xl_min, xl_max),
     legend.width=.35, legend.shrink=0.5, axis.args=list(cex.axis=.8), add = T)
dev.off()


