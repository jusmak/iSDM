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
scenario_labels = c('Interpolation', 'Extrapolation', 'Extrapolation2')

model_ind_all <- c(4:7,10)

#predictive accuracy
#interpolation
auc <- matrix(NA, length(species_list), length(model_ind_all))
tss <- matrix(NA, length(species_list), length(model_ind_all))
kappa <- matrix(NA, length(species_list), length(model_ind_all))
sens <- matrix(NA, length(species_list), length(model_ind_all))
spec <- matrix(NA, length(species_list), length(model_ind_all))
log_lik <- matrix(NA, length(species_list), length(model_ind_all))
tjur <- matrix(NA, length(species_list), length(model_ind_all))

for (s in 1:length(species_list)) {
  
  .species <- species_list[s]

  # Get validation results
  load(file = paste0('Model_fits_v', model_v, '/INLA_', .species, '_validation_', thin_mark, thin_PA_mark, '_quad_n_',
                    target_n_obs, '_', weights_mark, '_inla.RData'))
  
  #interpolation    
  scenario_ind <- grep('Int', names(output_validation_cv))

  for (j in 1:length(model_ind_all)) {
    j_temp <- model_ind_all[j]
    auc[s,j] <- round(mean(unlist(lapply(lapply(output_validation_cv, function(x)  unlist(x[j_temp])), function(x) x[[1]]))[scenario_ind]),2)
    kappa[s,j] <- round(mean(unlist(lapply(lapply(output_validation_cv, function(x)  unlist(x[j_temp])), function(x) x[[2]]))[scenario_ind]),2)
  }
}

model_label_temp <- rep(model_labels[model_ind_all[2]],length(species_list))
for (i in c(3:5)) {
  model_label_temp <- c(model_label_temp, rep(model_labels[model_ind_all[i]],length(species_list)))
}

model_order <- model_labels[model_ind_all[c(2:5)]]

#difference between models
auc_temp <- auc_temp <- auc[,c(2:5)] - auc[,1]
auc_int <- data.frame(cbind(as.vector(auc_temp), model_label_temp))
colnames(auc_int) <- c('AUC', 'Model')
auc_int$AUC <- as.numeric(auc_int$AUC)

#difference between models
kappa_temp <- kappa[,c(2:5)] - kappa[,1]
kappa_int <- data.frame(cbind(as.vector(kappa_temp), model_label_temp))
colnames(kappa_int) <- c('Kappa', 'Model')
kappa_int$Kappa <- as.numeric(kappa_int$Kappa)

#absolute model value
auc_abs <- data.frame(cbind(as.vector(auc[,1]), rep(model_labels[model_ind_all[1]],length(species_list))))
colnames(auc_abs) <- c('AUC', 'Model')
auc_abs$AUC <- as.numeric(auc_abs$AUC)

#absolute model value
kappa_abs <- data.frame(cbind(as.vector(kappa[,1]), rep(model_labels[model_ind_all[1]],length(species_list))))
colnames(kappa_abs) <- c('Kappa', 'Model')
kappa_abs$Kappa <- as.numeric(kappa_abs$Kappa)

g1 <- list()
g1[[1]] <- ggplot(data=auc_int, aes(x=factor(Model, level = model_order), y=AUC, fill=factor(Model, level = model_order))) + 
  geom_hline(yintercept = 0, size = 1.5, color = 'grey') +
  geom_boxplot(alpha = .75) +
  scale_fill_manual(values = c(wes_palette('Darjeeling2', n = 5)[2:4], wes_palette('Moonrise3', n = 5)[3])) +
  stat_summary(fun.data='mean_sdl', geom='crossbar', width=0.06, fill='white') +
  #ylim(-.25,.25) +
  coord_cartesian(ylim = c(-.1,.35)) +
  scale_y_continuous(breaks = round(seq(-.3 , .3, .1),1)) +
  #theme(axis.text.x = element_blank()) +
  ylab(expression(paste(Delta, ' AUC'))) +
  xlab('') +
  theme_light() +
  theme(legend.position = 'none', axis.text.x = element_text(angle = 15, vjust = .5, hjust=.5))

g1[[2]] <- ggplot(data=kappa_int, aes(x=factor(Model, level = model_order), y=Kappa, fill=factor(Model, level = model_order))) + 
  geom_hline(yintercept = 0, size = 1.5, color = 'grey') +
  geom_boxplot(alpha = .75) +
  scale_fill_manual(values = c(wes_palette('Darjeeling2', n = 5)[2:4], wes_palette('Moonrise3', n = 5)[3])) +
  stat_summary(fun.data='mean_sdl', geom='crossbar', width=0.06, fill='white') +
  coord_cartesian(ylim = c(-.15,.45)) +
  scale_y_continuous(breaks = round(seq(-.4, .4, .2),1)) +
  ylab(expression(paste(Delta, ' Kappa'))) +
  xlab('') +
  theme_light() +
  theme(axis.text.x = element_blank()) +
  theme(legend.position = 'none')

g1[[3]] <- ggplot(data=auc_abs, aes(x=Model, y=AUC, fill=Model)) + 
  #geom_hline(yintercept = median(auc_abs$AUC), size = 1.5, color = 'grey') +
  geom_boxplot(alpha = .75) +
  scale_fill_manual(values = wes_palette('Darjeeling2', n = 5)[1]) +
  stat_summary(fun.data='mean_sdl', geom='crossbar', width=0.06, fill='white') +
  coord_cartesian(ylim = mean(auc_abs$AUC) + c(-.3,.3)) +
  #scale_y_continuous(breaks = seq(.7, 1, .1)) +
  #theme(axis.text.x = element_blank()) +
  xlab('') +
  ylab('AUC') +
  theme_light() +
  theme(legend.position = 'none', axis.text.x = element_text(angle = 90, vjust = .5, hjust=.5))

g1[[4]] <- ggplot(data=kappa_abs, aes(x=Model, y=Kappa, fill=Model)) + 
  #geom_hline(yintercept = median(kappa_abs$Kappa), size = 1.5, color = 'grey') +
  geom_boxplot(alpha = .75) +
  scale_fill_manual(values = wes_palette('Darjeeling2', n = 5)[1]) +
  stat_summary(fun.data='mean_sdl', geom='crossbar', width=0.06, fill='white') +
  coord_cartesian(ylim = mean(kappa_abs$Kappa)+c(-.35,.35)) +
  #scale_y_continuous(breaks = seq(.1, .6, .1)) +
  xlab('') +
  ylab('Kappa') +
  theme_light() +
  theme(axis.text.x = element_blank()) +
  theme(legend.position = 'none')



#predictive accuracy
#extrapolation - 2
auc <- matrix(NA, length(species_list), length(model_ind_all))
tss <- matrix(NA, length(species_list), length(model_ind_all))
kappa <- matrix(NA, length(species_list), length(model_ind_all))
sens <- matrix(NA, length(species_list), length(model_ind_all))
spec <- matrix(NA, length(species_list), length(model_ind_all))
log_lik <- matrix(NA, length(species_list), length(model_ind_all))
tjur <- matrix(NA, length(species_list), length(model_ind_all))

for (s in 1:length(species_list)) {
  
  .species <- species_list[s]
  # Get validation results
  load(file = paste0('Model_fits_v', model_v, '/INLA_', .species, '_validation_', thin_mark, thin_PA_mark, '_quad_n_',
                     target_n_obs, '_', weights_mark, '_inla.RData'))
  
  #extrapolation    
  scenario_ind <- grep('Ext', names(output_validation_cv))
  
  for (j in 1:length(model_ind_all)) {
    j_temp <- model_ind_all[j]
    auc[s,j] <- round(mean(unlist(lapply(lapply(output_validation_cv, function(x)  unlist(x[j_temp])), function(x) x[[1]]))[scenario_ind]),2)
    kappa[s,j] <- round(mean(unlist(lapply(lapply(output_validation_cv, function(x)  unlist(x[j_temp])), function(x) x[[2]]))[scenario_ind]),2)
  }
}

#difference between models
auc_temp <- auc_temp <- auc[,c(2:5)] - auc[,1]
auc_int <- data.frame(cbind(as.vector(auc_temp), model_label_temp))
colnames(auc_int) <- c('AUC', 'Model')
auc_int$AUC <- as.numeric(auc_int$AUC)

#difference between models
kappa_temp <- kappa[,c(2:5)] - kappa[,1]
kappa_int <- data.frame(cbind(as.vector(kappa_temp), model_label_temp))
colnames(kappa_int) <- c('Kappa', 'Model')
kappa_int$Kappa <- as.numeric(kappa_int$Kappa)

#absolute model value
auc_abs <- data.frame(cbind(as.vector(auc[,1]), rep(model_labels[model_ind_all[1]],length(species_list))))
colnames(auc_abs) <- c('AUC', 'Model')
auc_abs$AUC <- as.numeric(auc_abs$AUC)

#absolute model value
kappa_abs <- data.frame(cbind(as.vector(kappa[,1]), rep(model_labels[model_ind_all[1]],length(species_list))))
colnames(kappa_abs) <- c('Kappa', 'Model')
kappa_abs$Kappa <- as.numeric(kappa_abs$Kappa)


g2 <- list()
g2[[1]] <- ggplot(data=auc_int, aes(x=factor(Model, level = model_order), y=AUC, fill=factor(Model, level = model_order))) + 
  #geom_hline(yintercept = 0, size = 1.5, color = 'grey') +
  geom_boxplot(alpha = .75) +
  #scale_fill_manual(values = c(wes_palette('Darjeeling2', n = 5)[2:4], wes_palette('Moonrise3', n = 5)[3])) +
  scale_fill_manual(values = rep(wes_palette('Darjeeling2', n = 5)[1], 4)) +
  #stat_summary(fun.data='mean_sdl', geom='crossbar', width=0.06, fill='white') +
  coord_cartesian(ylim = c(-.1,.25)) +
  scale_y_continuous(breaks = round(seq(-.3 , .3, .1),1)) +
  ylab(expression(paste(Delta, ' AUC'))) +
  xlab('') +
  theme_light() +
  theme(legend.position = 'none', axis.text.x = element_blank())

g2[[2]] <- ggplot(data=kappa_int, aes(x=factor(Model, level = model_order), y=Kappa, fill=factor(Model, level = model_order))) + 
  #geom_hline(yintercept = 0, size = 1.5, color = 'grey') +
  geom_boxplot(alpha = .75) +
  #scale_fill_manual(values = c(wes_palette('Darjeeling2', n = 5)[2:4], wes_palette('Moonrise3', n = 5)[3])) +
  scale_fill_manual(values = rep(wes_palette('Darjeeling2', n = 5)[1], 4)) +
  #stat_summary(fun.data='mean_sdl', geom='crossbar', width=0.06, fill='white') +
  coord_cartesian(ylim = c(-.2,.45)) +
  scale_y_continuous(breaks = round(seq(-.6, .6, .2),1)) +
  ylab(expression(paste(Delta, ' Kappa'))) +
  xlab('') +
  theme_light() +
  theme(axis.text.x = element_blank()) +
  theme(legend.position = 'none')


g2[[3]] <- ggplot(data=auc_abs, aes(x=Model, y=AUC, fill=Model)) + 
  #geom_hline(yintercept = median(auc_abs$AUC), size = 1.5, color = 'grey') +
  geom_boxplot(alpha = .75) +
  scale_fill_manual(values = wes_palette('Darjeeling2', n = 5)[1]) +
  #stat_summary(fun.data='mean_sdl', geom='crossbar', width=0.06, fill='white') +
  coord_cartesian(ylim = mean(auc_abs$AUC) + c(-.3,.3)) +
  scale_y_continuous(breaks = seq(.4, 1, .1)) +
  ylab('AUC') +
  xlab('') +
  theme_light() +
  theme(legend.position = 'none', axis.text.x = element_blank())

g2[[4]] <- ggplot(data=kappa_abs, aes(x=Model, y=Kappa, fill=Model)) + 
  #geom_hline(yintercept = median(kappa_abs$Kappa), size = 1.5, color = 'grey') +
  geom_boxplot(alpha = .75) +
  scale_fill_manual(values = wes_palette('Darjeeling2', n = 5)[1]) +
  #stat_summary(fun.data='mean_sdl', geom='crossbar', width=0.06, fill='white') +
  coord_cartesian(ylim = mean(kappa_abs$Kappa) + c(-.3,.3)) +
  #scale_y_continuous(breaks = seq(0, .9, .1), minor_breaks = seq(0 , .9, .1)) +
  ylab('Kappa') +
  xlab('') +
  theme_light() +
  theme(axis.text.x = element_blank()) +
  theme(legend.position = 'none')


##Accuracy of range prediction
load(file = 'Model_fits_v1/Summaries_cv_sp/Range_map_accuracy_50km.RData')

auc = matrix(unlist(lapply(range_accuracy, function(x) x[1,model_ind_all])), byrow = T, ncol = length(model_ind_all))
kappa = matrix(unlist(lapply(range_accuracy, function(x) x[2,model_ind_all])), byrow = T, ncol = length(model_ind_all))
tss = matrix(unlist(lapply(range_accuracy, function(x) x[3,model_ind_all])), byrow = T, ncol = length(model_ind_all))
spec = matrix(unlist(lapply(range_accuracy, function(x) x[4,model_ind_all])), byrow = T, ncol = length(model_ind_all))
sens = matrix(unlist(lapply(range_accuracy, function(x) x[5,model_ind_all])), byrow = T, ncol = length(model_ind_all))

#difference between models
auc_temp <- auc_temp <- auc[,c(2:5)] - auc[,1]
auc_int <- data.frame(cbind(as.vector(auc_temp), model_label_temp))
colnames(auc_int) <- c('AUC', 'Model')
auc_int$AUC <- as.numeric(auc_int$AUC)

#difference between models
kappa_temp <- kappa[,c(2:5)] - kappa[,1]
kappa_int <- data.frame(cbind(as.vector(kappa_temp), model_label_temp))
colnames(kappa_int) <- c('Kappa', 'Model')
kappa_int$Kappa <- as.numeric(kappa_int$Kappa)

#absolute model value
auc_abs <- data.frame(cbind(as.vector(auc[,1]), rep(model_labels[model_ind_all[1]],length(species_list))))
colnames(auc_abs) <- c('AUC', 'Model')
auc_abs$AUC <- as.numeric(auc_abs$AUC)

#absolute model value
kappa_abs <- data.frame(cbind(as.vector(kappa[,1]), rep(model_labels[model_ind_all[1]],length(species_list))))
colnames(kappa_abs) <- c('Kappa', 'Model')
kappa_abs$Kappa <- as.numeric(kappa_abs$Kappa)


g3 <- list()
g3[[1]] <- ggplot(data=auc_int, aes(x=factor(Model, level = model_order), y=AUC, fill=factor(Model, level = model_order))) + 
  #geom_hline(yintercept = 0, size = 1.5, color = 'grey') +
  geom_boxplot(alpha = .75) +
  #scale_fill_manual(values = c(wes_palette('Darjeeling2', n = 5)[2:4], wes_palette('Moonrise3', n = 5)[3])) +
  scale_fill_manual(values = rep(wes_palette('Darjeeling2', n = 5)[1], 4)) +
  #stat_summary(fun.data='mean_sdl', geom='crossbar', width=0.06, fill='white') +
  coord_cartesian(ylim = c(-.3,.2)) +
  #scale_y_continuous(breaks = round(seq(-.3 , .3, .1),1)) +
  ylab(expression(paste(Delta, ' AUC'))) +
  xlab('') +
  theme_light() +
  theme(axis.text.x = element_blank()) +
  theme(legend.position = 'none')

g3[[2]] <- ggplot(data=kappa_int, aes(x=factor(Model, level = model_order), y=Kappa, fill=factor(Model, level = model_order))) + 
  #geom_hline(yintercept = 0, size = 1.5, color = 'grey') +
  geom_boxplot(alpha = .75) +
  #scale_fill_manual(values = c(wes_palette('Darjeeling2', n = 5)[2:4], wes_palette('Moonrise3', n = 5)[3])) +
  scale_fill_manual(values = rep(wes_palette('Darjeeling2', n = 5)[1], 4)) +
  #stat_summary(fun.data='mean_sdl', geom='crossbar', width=0.06, fill='white') +
  #coord_cartesian(ylim = c(-.2,.65)) +
  #scale_y_continuous(breaks = round(seq(-.6, .6, .2),1)) +
  ylab(expression(paste(Delta, ' Kappa'))) +
  xlab('') +
  theme_light() +
  theme(axis.text.x = element_blank()) +
  theme(legend.position = 'none')


g3[[3]] <- ggplot(data=auc_abs, aes(x=Model, y=AUC, fill=Model)) + 
  #geom_hline(yintercept = median(auc_abs$AUC), size = 1.5, color = 'grey') +
  geom_boxplot(alpha = .75) +
  scale_fill_manual(values = wes_palette('Darjeeling2', n = 5)[1]) +
  #stat_summary(fun.data='mean_sdl', geom='crossbar', width=0.06, fill='white') +
  coord_cartesian(ylim = mean(auc_abs$AUC) + c(-.3,.1)) +
  #scale_y_continuous(breaks = seq(.4, 1, .1)) +
  ylab('AUC') +
  xlab('') +
  theme_light() +
  theme(axis.text.x = element_blank()) +
  theme(legend.position = 'none')

g3[[4]] <- ggplot(data=kappa_abs, aes(x=Model, y=Kappa, fill=Model)) + 
  #geom_hline(yintercept = median(kappa_abs$Kappa), size = 1.5, color = 'grey') +
  geom_boxplot(alpha = .75) +
  scale_fill_manual(values = wes_palette('Darjeeling2', n = 5)[1]) +
  #stat_summary(fun.data='mean_sdl', geom='crossbar', width=0.06, fill='white') +
  #coord_cartesian(ylim = mean(kappa_abs$Kappa) + c(-.3,.3)) +
  #scale_y_continuous(breaks = seq(0, .9, .1), minor_breaks = seq(0 , .9, .1)) +
  ylab('Kappa') +
  xlab('') +
  theme_light() +
  theme(axis.text.x = element_blank()) +
  theme(legend.position = 'none')


#Add predictive uncertainty and proportion of population inside the expected range
##load the uncertainty
load(file = 'Model_fits_v1/Summaries_cv_sp/Species_var_pred.RData')

model_labels <- c('PO', 'PA', 'Joint', 'PO', 'PA', 'PO+PA', 'PO+PA+samp.', 'PO+restr.', 'PA+restr.', 'PO+PA+restr.', 'PO+PA+samp.')
model_ind_all <- c(4:7,10)
model_order <- model_labels[model_ind_all]

var_ext = (sp_var[[5]][, model_ind_all] + sp_var[[6]][, model_ind_all] +
  sp_var[[7]][, model_ind_all] + sp_var[[8]][, model_ind_all])/4


model_label_temp <- rep(model_labels[model_ind_all[2]],length(species_list))
for (i in c(3:5)) {
  model_label_temp <- c(model_label_temp, rep(model_labels[model_ind_all[i]],length(species_list)))
}
model_order <- model_labels[model_ind_all[c(2:5)]]

var_part_rel <- var_ext[,2:5] - var_ext[,1]
var_part_rel <- data.frame(cbind(as.vector(var_part_rel), model_label_temp))
colnames(var_part_rel) <- c('Sd', 'Model')
var_part_rel$Sd <- as.numeric(var_part_rel$Sd)

var_part_abs <- data.frame(cbind(as.vector(var_ext[,1]), rep(model_labels[model_ind_all[1]],length(species_list))))
colnames(var_part_abs) <- c('Sd', 'Model')
var_part_abs$Sd <- as.numeric(var_part_abs$Sd)


g4 <- vector(mode='list', length=2)
g4[[1]] <- ggplot(data=var_part_abs, aes(x=Model, y=Sd, fill=Model)) +
  #geom_hline(yintercept = median(var_part_abs$Sd), size = 1.5, color = 'grey') +
  geom_boxplot(alpha = .75) +
  scale_fill_manual(values = wes_palette('Darjeeling2', n = 5)[1]) +
  #stat_summary(fun.data='mean_sdl', geom='crossbar', width=0.06, fill='white') +
  #coord_cartesian(ylim = mean(var_part_abs$Sd) + c(-4,4)) +
  ylab('Standard deviation') +
  xlab('') +
  theme_light() +
  theme(legend.position = 'none')
  

g4[[2]] <- ggplot(data=var_part_rel, aes(x=factor(Model, level = model_order), y=Sd, fill=factor(Model, level = model_order))) + 
  #geom_hline(yintercept = 0, size = 1.5, color = 'grey') +
  geom_boxplot(alpha = .75) +
  #scale_fill_manual(values = c(wes_palette('Darjeeling2', n = 5)[2:4], wes_palette('Moonrise3', n = 5)[3])) +
  scale_fill_manual(values = rep(wes_palette('Darjeeling2', n = 5)[1], 4)) +
  #stat_summary(fun.data='mean_sdl', geom='crossbar', width=0.06, fill='white') +
  #coord_cartesian(ylim = c(-2,1)) +
  ylab(expression(paste(Delta, ' Standard deviation'))) +
  xlab('') +
  theme_light() +
  theme(legend.position = 'none', axis.text.x = element_text(angle = 15, vjust = .5, hjust=.5))



#proportion of species densities falling inside the species range
model_labels <- c('PO', 'PA', 'PO+PA', 'PO', 'PA', 'PO+PA', 'PO+PA+samp.', 'PO+restr.', 'PA+restr.', 'PO+PA+restr.', 'PO+PA+samp.')
model_ind_all <- c(4:7,10)
model_order <- model_labels[model_ind_all]

load(file = 'Model_fits_v1/Summaries_cv_sp/Dens_range_map.RData')

model_label_temp <- rep(model_labels[model_ind_all[2]],length(species_list))
for (i in c(3:5)) {
  model_label_temp <- c(model_label_temp, rep(model_labels[model_ind_all[i]],length(species_list)))
}

model_order <- model_labels[model_ind_all[c(2:5)]]
dens_range$range = dens_range$range[, model_ind_all]
pop_rel <- dens_range$range[,2:5] - dens_range$range[,1]
pop_rel <- data.frame(cbind(as.vector(pop_rel), model_label_temp))
colnames(pop_rel) <- c('Pop', 'Model')
pop_rel$Pop <- as.numeric(pop_rel$Pop)

pop_abs <- data.frame(cbind(as.vector(dens_range$range[,1]), rep(model_labels[model_ind_all[1]],length(species_list))))
colnames(pop_abs) <- c('Pop', 'Model')
pop_abs$Pop <- as.numeric(pop_abs$Pop)

g5 <- vector(mode='list', length=2)
g5[[1]] <- ggplot(data=pop_abs, aes(x=Model, y=Pop, fill=Model)) +
  #geom_hline(yintercept = median(pop_abs$Pop), size = 1.5, color = 'grey') +
  geom_boxplot(alpha = .75) +
  scale_fill_manual(values = wes_palette('Darjeeling2', n = 5)[1]) +
  #stat_summary(fun.data='mean_sdl', geom='crossbar', width=0.06, fill='white') +
  coord_cartesian(ylim = mean(pop_abs$Pop) + c(-.6,.6)) +
  ylab('Proportion') +
  xlab('') +
  theme_light() +
  theme(legend.position = 'none')

g5[[2]] <- ggplot(data=pop_rel, aes(x=factor(Model, level = model_order), y=Pop, fill=factor(Model, level = model_order))) + 
  #geom_hline(yintercept = 0, size = 1.5, color = 'grey') +
  geom_boxplot(alpha = .75) +
  #scale_fill_manual(values = c(wes_palette('Darjeeling2', n = 5)[2:4], wes_palette('Moonrise3', n = 5)[3])) +
  scale_fill_manual(values = rep(wes_palette('Darjeeling2', n = 5)[1], 4)) +
  #stat_summary(fun.data='median_sdl', geom='crossbar', width=0.06, fill='white') +
  coord_cartesian(ylim = c(-1,.4)) +
  ylab(expression(paste(Delta, ' Proportion'))) +
  xlab('') +
  theme_light() +
  theme(legend.position = 'none', axis.text.x = element_text(angle = 15, vjust = .5, hjust=.5))


pdf('Figures/Fig2.pdf', width=8, height=7)
grid.arrange(g2[[3]], g2[[1]], grob(NULL), g2[[4]], g2[[2]],
             g3[[3]], g3[[1]], grob(NULL), g3[[4]], g3[[2]], 
             g5[[1]], g5[[2]], grob(NULL), g4[[1]], g4[[2]], nrow = 3, widths = c(1.5,3.5,.3,1.5,3.5))
dev.off()

