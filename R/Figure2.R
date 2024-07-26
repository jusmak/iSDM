## Manuscript: Fig 2

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


model_labels = c('PO', 'PA', 'PO+PA', 'PO+PA+samp.', 'PO+PA+restr.')
model_ind_all = c(4:7,10)

cov_labels = c('Temp', 'Temp_range', 'Prec', 'Prec_season', 'Cloud', 'EVI', 'TRI')
cov_labels = c(cov_labels, paste0(cov_labels, '_sqrd'))

cov_labels_plot = c('Temperature', 'Temperature range', 'Precipitation', 'Precipitation seasonality', 'Cloud coverage', 'EVI', 'TRI')

#get beta estimates
load(file = paste(.wd, "/Manuscript/Figures_manuscript_v12/Beta_estimates.RData", sep = ''))
beta_mean = beta_estimates$mean


#beta estimates
beta_mean = beta_mean[,,model_ind_all]

n_model = length(model_ind_all)
n_cov = length(cov_labels)
beta_array = array(beta_mean)
cov_array = rep(cov_labels, n_model, each = dim(beta_mean)[1])
model_array = rep(model_labels, each = n_cov*dim(beta_mean)[1])
beta_df = data.frame(beta_array, cov_array, model_array)
colnames(beta_df) = c('Beta', 'Covariate', 'Model')


#beta estimates for temperature, precipitation and EVI
g1 = list()

for (i in 1:(n_cov/2)) {
  cov_ind = c(i, i + (n_cov/2))
  
  beta_df_temp = beta_df[(beta_df$Covariate==cov_labels[i] | beta_df$Covariate==cov_labels[i + (n_cov/2)]),]
  
  g1[[i]] = ggplot() +
    geom_boxplot(data=beta_df_temp, aes(y=Beta, fill=Covariate, x=factor(Model, level = model_labels))) +
    theme_light() +
    ggtitle(cov_labels_plot[i]) +
    theme(legend.position = "none", text = element_text(size = 8), plot.title = element_text(size=8)) +
    # coord_cartesian(ylim = c(-12,12)) +
    scale_fill_manual(name = "Effect", labels = c("First order", "Second order"), values = wes_palette("Darjeeling2", n = 2)) +
    xlab('')
}

#remove y-axis label
for (i in 2:length(g1)) {
  g1[[i]] = g1[[i]] +
    ylab('')
}

#remove x-axis tick label
for (i in 1:(length(g1)-2)) {
  g1[[i]] = g1[[i]] +
    theme(axis.text.x = element_blank())
}


g1[[6]] = g1[[6]] +
  theme(legend.position = c(.5, -.3), legend.direction = 'horizontal', plot.margin = unit(c(5.5,5.5,25,5.5), 'points'))


in2mm = 25.4
pdf('Figures/Fig2.pdf', width = 168/in2mm, height = 160/in2mm)
grid.arrange(g1[[1]], g1[[3]], g1[[6]], nrow = 3, heights = c(1,1,1.2))
dev.off()

png('Figures/Fig2.png', width = 168/in2mm, height = 160/in2mm, units = 'in', res = 600)
grid.arrange(g1[[1]], g1[[3]], g1[[6]], nrow = 3, heights = c(1,1,1.2))
dev.off()
