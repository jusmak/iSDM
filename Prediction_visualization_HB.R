Prediction_visualization <- function(wd, model_fits, pred_data, n_samp, species, scale_out) {
  
  options(mc.cores = parallel::detectCores())
  
  #visualize the species observations and the offsets
  #load raster
  raster_cov_temp <- stack(paste(wd,"Data/Environment/Chelsa_SA.tif", sep = '/'))
  
  #load geographic extent
  domain_temp <- stack(paste(wd,"/Data/Domains/", species, ".tif", sep = ''))
  
  raster_cov_temp <- crop(raster_cov_temp, domain_temp)
  values(raster_cov_temp)[is.na(values(domain_temp)[,2]),] <- NA

  #rescale raster
  raster_cov_temp <- raster::aggregate(raster_cov_temp, fact = scale_out)

  #find cells which have all covariate values
  cov_df <- as.data.frame(raster_cov_temp, na.rm = FALSE)
  rows_na <- apply(cov_df,1,anyNA)
  ind_na <- which(rows_na==FALSE)

  #create a background
  background <- raster_cov_temp$Chelsa_SA.1
  sp_occ <- model_fits$training_data$coordinates[model_fits$training_data$response==1,]
  sp_occ <- cbind(sp_occ[,1]*1000 + model_fits$training_data$min_coordinates[1],
                  sp_occ[,2]*1000 + model_fits$training_data$min_coordinates[2])
  
  raster_cov_temp <- stack(paste(wd,"Data/Environment/Chelsa_SA.tif", sep = '/'))
  extent_temp <- c(-4e6, 5.5e6, -7.5e6, 2e6)
  raster_cov_temp <- crop(raster_cov_temp, extent_temp)
  background_domain <- raster_cov_temp$Chelsa_SA.1
  
  
  mycol <- rgb(0, 0, 255, max = 255, alpha = 75, names = "blue50")
  col_temp <- colorRamp("gray")
  mycol_2 <- col_temp(1)/255

  cm.cols1.1=function(x,bias=1) { colorRampPalette(c('grey90','grey90','grey90','grey90','grey60','grey60','grey60','steelblue4','steelblue1','gold','gold','red1','red1','red4','red4'),bias=bias)(x)
  }
  col_pal <- cm.cols1.1(50)
  
  #Figure 1: study domain and observations
  tiff(file = paste(wd, "/Figures/", species, "_Observations_offsets.tiff", sep = ''),
       width=10, height=10, units = "cm", res=400, type = 'cairo')
  
  par(mfrow = c(2,2), cex.main = .7, mai = c(0.1,0.1,0.2,0.1))
  values(background)[ind_na] <- 0
  values(background)[-ind_na] <- NA
  
  plot(background_domain, axes = F, box = FALSE, main = "Study domain", col = mycol_2, legend = FALSE)
  plot(background, axes = F, box = FALSE, legend = FALSE, add = TRUE, col = mycol)
  
  plot(background, axes = F, box = FALSE, legend = FALSE, main = paste("Observations, n_obs = ", as.character(sum(model_fits$training_data$response)), sep = ''),
       col = mycol)
  points(sp_occ, cex = .5, pch = 20)

  values(background)[ind_na] <- log(pred_data$offset[[2]][,1])
  plot(background, axes = F, box = FALSE, main = "Expert range map", xaxt='n', yaxt='n', col = col_pal,
       legend.args=list(text="log(lambda)", line=.2 ,side=3, cex=.5), axis.args=list(cex.axis=.5))

  values(background)[ind_na] <- log(pred_data$offset[[3]][,2])
  plot(background, axes = F, box = FALSE, main = "Elevation offset", xaxt='n', yaxt='n', col = col_pal,
       legend.args=list(text="log(lambda)", line=.2 ,side=3, cex=.5), axis.args=list(cex.axis=.5))

  dev.off()
  
  
  ## Hierarchical Bayesian ##
  
  #Figure 2: Predictions of each combination of offsets for fitting and predicting
  #take a random sample of draws
  set.seed(1)
  beta_samp <- extract(model_fits$HB_model[[1]])
  ind_samples <- sample(1:nrow(beta_samp$alpha), n_samp)

  #matrix of covariate values in prediction cells  
  pred_matrix <- cbind(pred_data$covariates, pred_data$covariates^2)
  
  tiff(file = paste(wd, "/Figures/", species, "_HB_mean_prediction.tiff", sep = ''),
       width=15, height=15, units = "cm", res=400, type = 'cairo')
  
  par(mfrow = c(4,4), cex.main = .7, mai = c(0.1,0.1,0.2,.1))
  title_temp <- c("none", "exp", "elev", "exp + elev")
  for (i in 1:length(model_fits)) {
    beta_samp <- extract(model_fits$HB_model[[i]])

    pred <- matrix(rep(beta_samp$alpha[ind_samples], nrow(pred_matrix)), nrow = nrow(pred_matrix), byrow = TRUE) + pred_matrix %*% t(beta_samp$beta[ind_samples,])
    for (j in 1:length(model_fits)) {
      offset_temp <- log(pred_data$offset[[j]][,1]) - mean(log(pred_data$offset[[j]][,1])) +
        log(pred_data$offset[[j]][,2]) - mean(log(pred_data$offset[[j]][,2]))
      mean_pred <- exp(pred + offset_temp)
      mean_pred <- apply(mean_pred, 1, mean)
      values(background)[ind_na] <- log(mean_pred)
      plot(background, axes = F, box = FALSE, legend = FALSE, main = title_temp[i], xaxt='n', yaxt='n', col = col_pal[1])
      plot(background, axes = F, box = FALSE, main = title_temp[i], xaxt='n', yaxt='n', col = col_pal, zlim = quantile(log(mean_pred), c(.2, 1)), add = TRUE,
           legend.args=list(text="log(lambda)", line=0 ,side=3, cex=.5), axis.args=list(cex.axis=.5),
           smallplot= c(.82,.84,.03,.2))
      
    }
  }
  
  dev.off()
  
  
  #Figure 3: Predictions where offsets are equal for fitting and predicting
  
  #compute thresholds
  source(paste(wd, "/R_code/Plot_thresholds_HB.R", sep = ''))
  thres_all <- Plot_thresholds(model_fits, ind_samples)
  thres_temp <- thres_all$TP05
  
  #Compile validation results
  pred_dens_cv <- matrix(NA, ncol = 2, nrow = length(model_fits$HB_model))
  auc_cv <-matrix(NA, ncol = 2, nrow = length(model_fits$HB_model))
  pred_dens_ind <- matrix(NA, ncol = 2, nrow = length(model_fits$HB_model))
  auc_ind <- matrix(NA, ncol = 2, nrow = length(model_fits$HB_model))
  pred_dens_train <- matrix(NA, ncol = 2, nrow = length(model_fits$HB_model))
  auc_train <- matrix(NA, ncol = 2, nrow = length(model_fits$HB_model))
  
  load(paste(.wd, '/', .species, '_validation.RData', sep = ''))
  
  for (i in 1:length(model_fits$HB_model)) {
    pred_dens_cv[i,1] <- pred_metric$HB_cv_valid$HB_cv_pred_dens[i,i]
    auc_cv[i,1] <- pred_metric$HB_cv_valid$HB_cv_auc[i,i]
    pred_dens_ind[i,1] <- pred_metric$HB_ind_valid$ind_validation_HB_pred_dens[i,i]
    auc_ind[i,1] <- pred_metric$HB_ind_valid$ind_validation_HB_auc[i,i]
    pred_dens_train[i,1] <- thres_all$pred_dens[i]
    auc_train[i,1] <- thres_all$AUC_int[i]    
  }
  validation_results <- list(auc_train, pred_dens_train,
                             auc_cv, pred_dens_cv, auc_ind,
                             pred_dens_ind)
  #rank models based on AUC
  for (i in seq(1,5,2)) {
    validation_results[[i]][,2] <- order(validation_results[[i]][,1], decreasing = TRUE)
    validation_results[[i]][,2] <- ifelse(validation_results[[i]][,2]==1,1,0)
  }
  #rank models based on predictive density
  for (i in seq(2,6,2)) {
    validation_results[[i]][,2] <- order(validation_results[[i]][,1], decreasing = FALSE)
    validation_results[[i]][,2] <- ifelse(validation_results[[i]][,2]==1,1,0)
  }
  
  #locations for validation figures
  y_max <- background@extent[4]
  y_int <- (background@extent[4]-background@extent[3])/background@nrows
  ats <- y_max - y_int * seq(1,6,1)*20
  

  pred_matrix <- cbind(pred_data$covariates, pred_data$covariates^2)
  
  tiff(file = paste(wd, "/Figures/", species, "_HB_mean_prediction_simple.tiff", sep = ''),
       width=17, height=5, units = "cm", res=400, type = 'cairo')
  
  par(mfrow = c(1,4), cex.main = .7, mai = c(0.1,0.1,0.2,.1))
  title_temp <- c("none", "exp", "elev", "exp + elev")
  
  for (i in 1:length(model_fits)) {
    beta_samp <- extract(model_fits$HB_model[[i]])
    
    pred <- matrix(rep(beta_samp$alpha[ind_samples], nrow(pred_matrix)), nrow = nrow(pred_matrix), byrow = TRUE) + pred_matrix %*% t(beta_samp$beta[ind_samples,])
    for (j in i) {
      offset_temp <- log(pred_data$offset[[j]][,1]) - mean(log(pred_data$offset[[j]][,1])) +
        log(pred_data$offset[[j]][,2]) - mean(log(pred_data$offset[[j]][,2]))
      mean_pred <- exp(pred + offset_temp)
      mean_pred <- apply(mean_pred, 1, mean)
      values(background)[ind_na] <- log(mean_pred)
      l_pred <- values(background)
      
      #set color scales
      in.max=max(l_pred,na.rm=T)
      in.min=quantile(l_pred,prob=1e-2,na.rm=T)
      
      tmp=quantile(l_pred[l_pred>=thres_temp[i]],seq(0,1,length=50),na.rm=T)
      hb.breaks=c(seq(in.min,thres_temp[i],length=25),tmp)
      values(background)[values(background)<in.min] <- in.min
      mean.col.brk <- hb.breaks
      
      col.brk=data.frame(cols=c(grey(seq(.8,.2,length=24)),
                                colorRampPalette(c("steelblue4", "steelblue1","gold","gold3", "red1", "red4"))(50),NA),
                         breaks=mean.col.brk,stringsAsFactors=F)
      col.brk=col.brk[!duplicated(col.brk$breaks),]
      
      #make sure the top color is red
      col.brk$cols[(nrow(col.brk)-1)]='#8B0000'      
      cex.factor <- 1
      t.ex <- ifelse(is.null(cex.factor),1,cex.factor)*.7
      scaleBarHeight <- .2
      
      #avoid too many legend labels
      lab.breaks.tmp=col.brk$breaks
      keep=round(seq(2,length(lab.breaks.tmp)-1,length=3))
      axis.labs=lab.breaks.tmp[keep]

      
      plot(background, axes = F, box = FALSE, main = title_temp[i], xaxt='n', yaxt='n',
           zlim = c(col.brk$breaks[1],col.brk$breaks[nrow(col.brk)]), col=col.brk$cols[-nrow(col.brk)],
           legend.args=list(text="log(lambda)", line=.2,side=3,cex=.7*t.ex),
           axis.args=list(at=axis.labs, labels=round(axis.labs,2),cex.axis=.8*t.ex),
           breaks=col.brk$breaks, smallplot= c(.72,.74,.03,scaleBarHeight))
      
      #add validation results
      mtext(substitute(paste('AUC'['train'],' = ',x1),
                       list(x1=round(validation_results[[1]][i,1],3))),
            cex=.3, 4,.5,las=1,adj=c(0,0), at = ats[1],
            col = ifelse(validation_results[[1]][i,2]==1, 'grey10','grey10'))
      mtext(substitute(paste('PD'['train'],' = ',x1),
                       list(x1=round(validation_results[[2]][i,1],1))),
            cex=.3, 4,.5,las=1,adj=c(0,0), at = ats[2],
            col = ifelse(validation_results[[2]][i,2]==1, 'grey10','grey10'))
      mtext(substitute(paste('AUC'['cv'],' = ',x1),
                       list(x1=round(validation_results[[3]][i,1],3))),
            cex=.3, 4,.5,las=1,adj = c(0,0), at = ats[3],
            col = ifelse(validation_results[[3]][i,2]==1, 'grey10','grey10'))
      mtext(substitute(paste('PD'['cv'],' = ',x1),
                       list(x1=round(validation_results[[4]][i,1],1))),
            cex=.3, 4,.5,las=1,adj = c(0,0), at = ats[4],
            col = ifelse(validation_results[[4]][i,2]==1, 'grey10','grey10'))
      mtext(substitute(paste('AUC'['ind.'],' = ',x1),
                       list(x1=round(validation_results[[5]][i,1],3))),
            cex=.3, 4,.5,las=1,adj=c(0,0), at = ats[5],
            col = ifelse(validation_results[[5]][i,2]==1, 'grey10','grey10'))
      mtext(substitute(paste('PD'['ind.'],' = ',x1),
                       list(x1=round(validation_results[[6]][i,1],1))),
            cex=.3, 4,.5,las=1,adj=c(0,0), at = ats[6],
            col = ifelse(validation_results[[6]][i,2]==1, 'grey10','grey10'))
      
    }
  }
  
  dev.off()
  
  #Figure 4: standard deviation of the predictive posterior distribution
  
  tiff(file = paste(wd, "/Figures/", species, "_HB_sd_prediction.tiff", sep = ''),
       width=17, height=5, units = "cm", res=400, type = 'cairo')
  
  par(mfrow = c(1,4), cex.main = .7, mai = c(0.1,0.1,0.2,.1))
  title_temp <- c("none", "exp", "elev", "exp + elev")
  for (i in 1:length(model_fits)) {
    beta_samp <- extract(model_fits$HB_model[[i]])
    pred <- matrix(rep(beta_samp$alpha, nrow(pred_matrix)), nrow = nrow(pred_matrix), byrow = TRUE) + pred_matrix %*% t(beta_samp$beta)
    for (j in i) {
      offset_temp <- log(pred_data$offset[[j]][,1]) - mean(log(pred_data$offset[[j]][,1])) +
        log(pred_data$offset[[j]][,2]) - mean(log(pred_data$offset[[j]][,2]))
      mean_pred <- pred + offset_temp
      mean_pred <- apply(mean_pred, 1, function(x) sd(x))
      values(background)[ind_na] <- mean_pred
      l_pred <- values(background)
      
      #set color scales
      in.max=quantile(l_pred,prob = .99, na.rm=T)
      in.min=min(l_pred,na.rm=T)
      
      thres_temp <- quantile(l_pred, prob = .25, na.rm=T)
      
      tmp=quantile(l_pred[l_pred>=thres_temp],seq(0,1,length=50),na.rm=T)
      hb.breaks=c(seq(in.min,thres_temp,length=25),tmp)
      values(background)[values(background)>in.max] <- in.max
      mean.col.brk <- hb.breaks
      
      col.brk=data.frame(cols=c(grey(seq(.8,.2,length=24)),
                                colorRampPalette(c("steelblue4", "steelblue1","gold","gold3", "red1", "red4"))(50),NA),
                         breaks=mean.col.brk,stringsAsFactors=F)
      col.brk=col.brk[!duplicated(col.brk$breaks),]
      
      #make sure the top color is red
      col.brk$cols[(nrow(col.brk)-1)]='#8B0000'
      
      #avoid too many legend labels
      lab.breaks.tmp=col.brk$breaks
      keep=round(seq(2,length(lab.breaks.tmp)-1,length=3))
      axis.labs=lab.breaks.tmp[keep]
      
      plot(background, axes = F, box = FALSE, main = title_temp[i], xaxt='n', yaxt='n',
           zlim = c(col.brk$breaks[1],col.brk$breaks[nrow(col.brk)]), col=col.brk$cols[-nrow(col.brk)],
           legend.args=list(text="sd", line=.2,side=3,cex=.7*t.ex),
           axis.args=list(at=axis.labs, labels=round(axis.labs,2),cex.axis=.8*t.ex),
           breaks=col.brk$breaks, smallplot= c(.72,.74,.03,scaleBarHeight))
    }
  }
  dev.off()
   
}