## Manuscript: Fig 3

#---- Load packages ----#
suppressPackageStartupMessages({
  library(raster)
  library(ggplot2)
  library(dplyr)
  library(pROC)
  library(INLA)
  library(gridExtra)
  library(gtable)
  library(RColorBrewer)
  library(wesanderson)
  library(tidyverse)
  library(maptools)
  library(rgdal)
  library(pracma)
  library(rnaturalearth)
  library(sf)
})



#---- Set working directory ----#
setwd('iSDM')
.wd = getwd()

#---- Initialize Environment ----#
.seed = 5326
set.seed(.seed)
options(mc.cores = parallel::detectCores())


model_labels = c('PO', 'PA', 'PO+PA', 'PO+PA+samp.', 'PO+PA+restr.')
model_ind_all = c(4:7,10)


#get predicted range map size
load(file = paste(.wd, "/Manuscript/Figures_manuscript_v12/All_sp/Species_range_pred_50km.RData", sep = ''))
load(file = paste(.wd, "/Manuscript/Figures_manuscript_v12/All_sp/Species_richness_range_maps_50km.RData", sep = ''))
load(file = paste(.wd, "/Manuscript/Figures_manuscript_v12/All_sp/Species_range_data.RData", sep = ''))


##Species richness prediction
#change the ocean cells to NA
env_raster = stack(paste(.wd,"Study_area/Study_area_laea.tif", sep = '/'))
#aggregate to 50km resolution
env_raster_agg = raster::aggregate(env_raster, fact = 50, fun = min)

#find cells which have NA value
ind_na = which(is.na(raster::values(env_raster_agg)))

#for plotting, change NAs to zeros
raster::values(env_raster_agg) = ifelse(!is.na(raster::values(env_raster_agg)), 1, 0)

#country borders
world = ne_countries(scale = "medium", returnclass = "sf")
world = st_transform(world, crs = crs(env_raster_agg))

#set the color scale for absolute values
col.brk=data.frame(cols=c(grey(seq(.8,.2,length=24)),colorRampPalette(c("steelblue4", "steelblue1","gold","gold3", "red1", "red4"))(50)))
#make sure the top color is red
col.brk$cols[(nrow(col.brk)-1)]='#8B0000'
#set the color scale for relative values
col.brk.rel=data.frame(cols=c(colorRampPalette(c("steelblue4", "steelblue1"))(15), grey(.9), colorRampPalette(c("gold3", "red1"))(15)))

#Coordinate labels
x = c(-85,-85,-68)
y = c(-7.5,0,10)
d = data.frame(lon=x, lat=y)
coordinates(d) = c("lon", "lat")
proj4string(d) = CRS("+init=epsg:4326") # WGS 84
CRS.new = raster::crs(env_raster_agg)
coord_trans = spTransform(d, CRS.new)
coord_temp = coord_trans@coords

g1 = list()

#map limits
xl_min = -1.3e6
xl_max = 1e6
yl_min = -.5e6
yl_max = 2e6

##Fig 1: Richness from expert range maps
temp_1 = species_rich
temp_1[ind_na] = NA

#transform to a data frame
test_spdf = as(temp_1, "SpatialPixelsDataFrame")
test_df = as.data.frame(test_spdf)
colnames(test_df) = c("value", "x", "y")

g1[[1]] = ggplot() +  
  geom_tile(data=test_df, aes(x=x, y=y, fill=value)) +
  geom_sf(data = world, fill = NA, size = .3, color = 'black') +
  coord_equal() +
  coord_sf(xlim = c(xl_min, xl_max), ylim = c(yl_min, yl_max)) +
  theme_light() +
  scale_fill_gradientn(colours = col.brk$cols) +
  scale_x_continuous(breaks=c(-85,-75,-67), labels=c('85\u00B0 W', '75\u00B0 W', '67\u00B0 W'),
                     minor_breaks = NULL) +
  scale_y_continuous(breaks=c(-10,0,10), labels=c('10\u00B0 S', '0\u00B0 N', '10\u00B0 N'),
                     minor_breaks = NULL) +
  ylab('Lat') +
  xlab('Lon') +
  ggtitle('Expert range map') +
  theme(legend.position = 'none', panel.background = element_rect(fill = "lightblue"),
        text = element_text(size = 8), plot.title = element_text(size=8),
        plot.margin = unit(c(5.5,5.5,1,5.5), 'points'))

#model predictions: PO, PA
for (i in c(1,2)) {
  temp_1 = pred_range$richness_map[[model_ind_all[i]]]
  temp_1[ind_na] = NA
  
  #transform to a data frame
  test_spdf = as(temp_1, "SpatialPixelsDataFrame")
  test_df = as.data.frame(test_spdf)
  colnames(test_df) = c("value", "x", "y")
  
  g1[[i+1]] = ggplot() +  
    geom_tile(data=test_df, aes(x=x, y=y, fill=value)) +
    geom_sf(data = world, fill = NA, size = .3, color = 'black') +
    coord_equal() +
    coord_sf(xlim = c(xl_min, xl_max), ylim = c(yl_min, yl_max)) +
    theme_light() +
    scale_fill_gradientn(colours = col.brk$cols) +
    scale_x_continuous(breaks=c(-85,-75,-67), labels=NULL,
                       minor_breaks = NULL) +
    scale_y_continuous(breaks=c(-10,0,10), labels=NULL,
                       minor_breaks = NULL) +
    ylab('Lat') +
    xlab('Lon') +
    ggtitle(model_labels[i]) +
    theme(legend.position = 'none', panel.background = element_rect(fill = "lightblue"),
          text = element_text(size = 8), plot.title = element_text(size=8),
          plot.margin = unit(c(5.5,5.5,10,5.5), 'points'),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank())
  
}

#model predictions: PO+PA+samp.
i=4
temp_1 = pred_range$richness_map[[model_ind_all[i]]]
temp_1[ind_na] = NA

#transform to a data frame
test_spdf = as(temp_1, "SpatialPixelsDataFrame")
test_df = as.data.frame(test_spdf)
colnames(test_df) = c("value", "x", "y")

lim_z = c(0, max(test_df$value, na.rm = T))

g1[[4]] = ggplot() +  
  geom_tile(data=test_df, aes(x=x, y=y, fill=value)) +
  geom_sf(data = world, fill = NA, size = .3, color = 'black') +
  coord_equal() +
  coord_sf(xlim = c(xl_min, xl_max), ylim = c(yl_min, yl_max)) +
  theme_light() +
  scale_fill_gradientn(name = NULL, colours = col.brk$cols, breaks = lim_z, labels = c(0, 'max')) +
  scale_x_continuous(breaks=c(-85,-75,-67), labels=NULL,
                     minor_breaks = NULL) +
  scale_y_continuous(breaks=c(-10,0,10), labels=NULL,
                     minor_breaks = NULL) +
  ylab('Lat') +
  xlab('Lon') +
  ggtitle(model_labels[i]) +
  theme(panel.background = element_rect(fill = "lightblue"),
        text = element_text(size = 8), plot.title = element_text(size=8),
        plot.margin = unit(c(5.5,5.5,10,5.5), 'points'),
        legend.text = element_text(size=8),
        legend.direction = 'horizontal',
        legend.position = c(.5,-.15),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        legend.key.height= unit(.05, 'in'),
        legend.key.width= unit(.15, 'in'))


#Observation density
load(file = paste(.wd, "/Manuscript/Figures_manuscript_v12/All_sp/Species_po_data.RData", sep = ''))

temp_obs = pred_range$richness_map[[model_ind_all[1]]]
temp_obs[ind_na] = NA
#raster::values(temp_obs)[-ind_na] = 0
#temp_obs = aggregate(temp_obs, fact=50)
counts = log(table(raster::cellFromXY(temp_obs,PO_temp)))
#counts[counts>quantile(counts,.995)] = quantile(counts,.995)
values(temp_obs)[as.numeric(names(counts))] = counts
values(temp_obs)[-as.numeric(names(counts))] = 0
temp_obs[ind_na] = NA

test_spdf = as(temp_obs, "SpatialPixelsDataFrame")
test_df = as.data.frame(test_spdf)
colnames(test_df) = c("value", "x", "y")

lim_z = c(0, max(test_df$value, na.rm = T))

g1[[5]] = ggplot() +  
  geom_tile(data=test_df, aes(x=x, y=y, fill=value)) +
  geom_sf(data = world, fill = NA, size = .3, color = 'black') +
  coord_equal() +
  coord_sf(xlim = c(xl_min, xl_max), ylim = c(yl_min, yl_max)) +
  theme_light() +
  scale_fill_gradientn(name = NULL, colours = brewer.pal(n = 9, name = "YlGnBu")[3:9]) +
  scale_x_continuous(breaks=c(-85,-75,-67), labels=NULL,
                     minor_breaks = NULL) +
  scale_y_continuous(breaks=c(-10,0,10), labels=NULL,
                     minor_breaks = NULL) +
  ylab('Lat') +
  xlab('Lon') +
  ggtitle('Sampling intensity') +
  theme(panel.background = element_rect(fill = "lightblue"),
        text = element_text(size = 8), plot.title = element_text(size=8),
        plot.margin = unit(c(5.5,5.5,10,5.5), 'points'),
        legend.text = element_text(size=8),
        legend.direction = 'horizontal',
        legend.position = c(.5,-.15),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        legend.key.height= unit(.05, 'in'),
        legend.key.width= unit(.15, 'in'))


in2mm = 25.4
pdf('Figures/Fig3.pdf', width = 168/in2mm, height = 45/in2mm)
grid.arrange(g1[[1]], g1[[2]], g1[[3]], g1[[4]], g1[[5]], nrow = 1, widths = c(1.25,1,1,1,1))
dev.off()

png('Figures/Fig3.png', width = 168/in2mm, height = 45/in2mm, units = 'in', res = 1200)
grid.arrange(g1[[1]], g1[[2]], g1[[3]], g1[[4]], g1[[5]], nrow = 1, widths = c(1.25,1,1,1,1))
dev.off()
