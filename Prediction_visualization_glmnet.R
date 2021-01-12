Prediction_visualization <- function(wd, model_fits, training_data, offset_full, pred_data, n_samp, species, scale_out, pred_metric) {
  
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
  values(background)[ind_na] <- 0
  values(background)[-ind_na] <- NA
  
  #define colors
  mycol <- rgb(0, 0, 255, max = 255, alpha = 75, names = "blue50")
  col_temp <- colorRamp("gray")
  mycol_2 <- col_temp(1)/255

  cm.cols1.1=function(x,bias=1) { colorRampPalette(c('grey90','grey90','grey90','grey90','grey60','grey60','grey60','steelblue4','steelblue1','gold','gold','red1','red1','red4','red4'),bias=bias)(x)
  }
  col_pal <- cm.cols1.1(50)
  
  ## GLMNET ##
  
  #Figure 2: Predictions of each combination of offsets for fitting and predicting
  
  tiff(file = paste(wd, "/Figures/", species, "_glmnet_mean_prediction.tiff", sep = ''),
       width=15, height=15, units = "cm", res=400, type = 'cairo')
  
  par(mfrow = c(4,4), cex.main = .7, mai = c(0.1,0.1,0.2,.1))
  title_temp <- c("none", "exp", "elev", "exp + elev")
  x_temp <- model.matrix(~ poly(pred_data$covariates, degree = 2, raw = TRUE))
  for (i in 1:length(model_fits$glmnet_model)) {
    glmnet_temp <- model_fits$glmnet_model[[i]]
    for (j in 1:length(model_fits$glmnet_model)) {
      offset_temp <- log(pred_data$offset[[j]][,1]) + log(pred_data$offset[[j]][,2])
      pred <- predict(object = glmnet_temp, newx = x_temp,
                      newoffset = offset_temp)
      
      mean_pred <- pred[,ncol(pred)]
      values(background)[ind_na] <- mean_pred
      plot(background, axes = F, box = FALSE, legend = FALSE, main = title_temp[i], xaxt='n', yaxt='n', col = col_pal[1])
      plot(background, axes = F, box = FALSE, main = title_temp[i], xaxt='n', yaxt='n', col = col_pal, zlim = quantile(mean_pred, c(.2, 1)), add = TRUE,
           legend.args=list(text="log(lambda)", line=0 ,side=3, cex=.5), axis.args=list(cex.axis=.5),
           smallplot= c(.82,.84,.03,.2))
    }
  }
  
  dev.off()
  
  #Figure 3: Predictions where offsets are equal for fitting and predicting
  
  #compute thresholds
  source(paste(wd, "/R_code/Plot_thresholds_glmnet.R", sep = ''))
  thres_all <- Plot_thresholds(model_fits, training_data, offset_full)
  thres_temp <- thres_all$TP05

  #Compile validation results
  pred_dens_train <- matrix(NA, ncol = 2, nrow = length(model_fits$HB_model))
  auc_train <- matrix(NA, ncol = 2, nrow = length(model_fits$HB_model))
  pred_dens_cv <- matrix(NA, ncol = 2, nrow = length(model_fits$HB_model))
  auc_cv <-matrix(NA, ncol = 2, nrow = length(model_fits$HB_model))
  pred_dens_ind <- matrix(NA, ncol = 2, nrow = length(model_fits$HB_model))
  auc_ind <- matrix(NA, ncol = 2, nrow = length(model_fits$HB_model))
  
  for (i in 1:length(model_fits$glmnet_model)) {
    pred_dens_train[i,1] <- pred_metric$glmnet_train_valid$train_validation_glmnet_pred_dens[i,i]
    auc_train[i,1] <- pred_metric$glmnet_train_valid$train_validation_glmnet_auc[i,i]
    pred_dens_cv[i,1] <- pred_metric$glmnet_cv_valid$glmnet_cv_pred_dens[i,i]
    auc_cv[i,1] <- pred_metric$glmnet_cv_valid$glmnet_cv_auc[i,i]
    pred_dens_ind[i,1] <- pred_metric$glmnet_ind_valid$ind_validation_glmnet_pred_dens[i,i]
    auc_ind[i,1] <- pred_metric$glmnet_ind_valid$ind_validation_glmnet_auc[i,i]
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
  
  tiff(file = paste(wd, "/Figures/", species, "_glmnet_mean_prediction_simple.tiff", sep = ''),
       width=17, height=5, units = "cm", res=400, type = 'cairo')
  
  par(mfrow = c(1,4), cex.main = .7, mai = c(0.1,0.1,0.2,.1))
  title_temp <- c("none", "exp", "elev", "exp + elev")
  
  x_temp <- model.matrix(~ poly(pred_data$covariates, degree = 2, raw = TRUE))
  
  for (i in 1:length(model_fits$glmnet_model)) {
    glmnet_temp <- model_fits$glmnet_model[[i]]
    for (j in i) {
      offset_temp <- log(pred_data$offset[[j]][,1]) + log(pred_data$offset[[j]][,2])
      pred <- predict(object = glmnet_temp, newx = x_temp,
                      newoffset = offset_temp)
      mean_pred <- pred[,ncol(pred)]
      values(background)[ind_na] <- mean_pred
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
  
}