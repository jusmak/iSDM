Inference_ppml_cv_x_fold <- function(wd, training_data, offset_full, ind_fold) {
  
  # output matrices
  ppml_cv_pred_dens <- matrix(NA, ncol = length(offset_full), nrow = length(offset_full))
  ppml_cv_auc <- matrix(NA, ncol = length(offset_full), nrow = length(offset_full))
  model_fits_cv <- list()
  
  # loop over different offsets used for model fitting
  for (i in 1:length(offset_full)) {
    
    # temporary output list and matrices for the specific offset
    # columns are for the offsets used for predicting and rows for cross-validation folds
    model_fits_cv_offset <- list()
    lik_folds_pred_dens <- matrix(NA, nrow = length(ind_fold), ncol = length(offset_full))
    lik_folds_auc <- matrix(NA, nrow = length(ind_fold), ncol = length(offset_full))
    
    for (j in 1:length(ind_fold)) {
      
      # keep other occurrence points as training data
      ind_temp <- ind_fold[[j]]
      ind_temp <- (1:sum(training_data$response==1))[-ind_temp]
      
      # keep all quadrature points in the training data
      ind_keep_1 <- c(ind_temp, (sum(training_data$response==1)+1):length(training_data$response))
      
      ##ADJUST DATA FOR THE SUBSAMPLE##
      
      #scale covariates to their original scales
      cov_temp <- data.frame(training_data$covariates)
      cov_temp_orig <- data.frame(mapply(function(x,y,z) (x*y+z), cov_temp, training_data$cov_sd, training_data$cov_mean))
      
      #make a subsample of the data sets
      cov_temp <- cov_temp_orig[ind_keep_1,]
      
      #scale covariates to have mean zero and sd one
      cov_mean_temp <- apply(cov_temp, 2, mean)
      cov_sd_temp <- apply(cov_temp, 2, sd)
      cov_temp_sc <- mapply(function(x,y,z) ((x-y)/z), cov_temp, cov_mean_temp, cov_sd_temp)
      
      #take a subsample of the occurrence records
      resp_temp <- training_data$response[ind_keep_1]
      
      #take a subsample of the weights
      weights_temp <- training_data$weights[ind_keep_1]
      
      #take a subsample of coordinates
      coord_temp <- training_data$coordinates[ind_keep_1,]
      
      #fine scale quadrature weights which are not related to occurrence record have weight one
      ind_overlap <- matrix(NA,  sum(resp_temp))
      ind_quad <- which(weights_temp == .5 & resp_temp == 0)
      for (l in 1:sum(resp_temp)) {
        ind_overlap[l] <- which.min((coord_temp[ind_quad,1]-coord_temp[l,1])^2+(coord_temp[ind_quad,2]-coord_temp[l,2])^2)
      }
      weights_temp[ind_quad[-ind_overlap]] <- 1
      
      #scale offsets to the number of occurrence points
      #scale offsets to their original scale
      offset_temp <- offset_full[[i]][ind_keep_1,]*(sum(training_data$response)/sum(training_data$weights))
      
      #compute population in the cross-validation fold
      pop_temp <- matrix(apply(offset_temp, 2, function(x) sum(x*weights_temp)), nrow = 1)
      
      #scale intensities to sum the population
      offset_temp <- mapply(function(x,y) x*sum(resp_temp)/y, data.frame(offset_temp), data.frame(pop_temp))
      
      #scale intensities to diverge from constant distribution
      offset_temp_1 <- offset_temp/(sum(resp_temp)/sum(weights_temp))

      form = ~ poly(cov_temp_sc, degree = 2, raw = TRUE)
      
      #set data as a data frame in the ppmlasso format
      #coordinates, covariates, response and weights
      #weights are products of the areas and the offsets
      ppm_data <- data.frame(cbind(coord_temp, cov_temp_sc, resp_temp,
                                   weights_temp * offset_temp_1[,1] * offset_temp_1[,2]))
      colnames(ppm_data) <- c('X', 'Y', colnames(training_data$covariates), 'Pres', 'wt')
      
      fit_ppml = ppmlasso(form, data = ppm_data, n.fits = 50, max.it = 10)
      #model feature is too large to be stored
      #model_fits_cv_offset[[j]] <- fit_ppml

      # index for prediction
      ind_temp <- ind_fold[[j]]
      ind_keep_2 <- c(ind_temp, (sum(training_data$response==1)+1):length(training_data$response))
      
      ##ADJUST DATA FOR THE SUBSAMPLE##
      
      #make a subsample of the data sets
      cov_temp_test <- cov_temp_orig[ind_keep_2,]
      
      #scale covariates to the training data scale
      cov_temp_test_sc <- mapply(function(x,y,z) ((x-y)/z), cov_temp_test, cov_mean_temp, cov_sd_temp)
      
      #take a subsample of the occurrence records
      resp_temp_test <- training_data$response[ind_keep_2]
      
      #take a subsample of the weights
      weights_temp_test <- training_data$weights[ind_keep_2]
      
      #take a subsample of coordinates
      coord_temp_test <- training_data$coordinates[ind_keep_2,]
      
      #fine scale quadrature weights which are not related to occurrence record have weight one
      ind_overlap <- matrix(NA,  sum(resp_temp_test))
      ind_quad <- which(weights_temp_test == .5 & resp_temp_test == 0)
      for (l in 1:sum(resp_temp_test)) {
        ind_overlap[l] <- which.min((coord_temp_test[ind_quad,1]-coord_temp_test[l,1])^2+(coord_temp_test[ind_quad,2]-coord_temp_test[l,2])^2)
      }
      weights_temp_test[ind_quad[-ind_overlap]] <- 1
      
      #loop over different offset choices for predictions
      for (k in 1:length(offset_full)) {
        
        #scale offsets to their original scale
        offset_temp <- offset_full[[k]][ind_keep_1,]*(sum(training_data$response)/sum(training_data$weights))
        
        #compute population in the cross-validation fold
        pop_temp <- matrix(apply(offset_temp, 2, function(x) sum(x*weights_temp)), nrow = 1)
        
        #scale intensities to sum the population
        offset_temp <- mapply(function(x,y) x*sum(resp_temp)/y, data.frame(offset_temp), data.frame(pop_temp))
        
        #get offset of the closest training point
        #offset for the quadrature points
        offset_temp_test <- offset_temp[resp_temp==0,]
        
        #offset for the presence points
        ind_close <- matrix(NA, nrow = sum(resp_temp_test==1))
        for (m in 1:length(ind_close)) {
          ind_close[m] <- which.min((coord_temp_test[m,1] - coord_temp[,1])^2 +
                                      (coord_temp_test[m,2] - coord_temp[,2])^2)
        }
        
        #compile offsets
        offset_temp_test <- rbind(offset_temp[ind_close,], offset_temp_test)
        offset_temp_test <- offset_temp_test/(sum(resp_temp)/sum(weights_temp))
        
        #compile a prediction data set
        ppm_data_pred <- data.frame(coord_temp_test, cov_temp_test_sc,
                          weights_temp_test * offset_temp_test[,1] * offset_temp_test[,2])
        
        colnames(ppm_data_pred) <- c('X', 'Y', colnames(training_data$covariates), 'wt')
        
        #predict
        pred_validation <- predict.ppmlasso(object = fit_ppml, newdata = ppm_data_pred)

        #predictive density
        pred_lik <- dpois(resp_temp_test, exp(pred_validation[,ncol(pred_validation)] + log(weights_temp_test)), log = TRUE)
        lik_folds_pred_dens[j,k] <- sum(pred_lik)
        
        #occurrence probability -> AUC
        occ_prob <- 1-exp(-exp(pred_validation[,ncol(pred_validation)] + log(weights_temp_test)))
        suppressMessages(
          myRoc <- pROC::roc(response = resp_temp_test, predictor = occ_prob))
        lik_folds_auc[j,k] <- myRoc$auc
      }
    }
    #store each model fit
    model_fits_cv[[i]] <- model_fits_cv_offset
    
    #compute average predictive density and AUC
    ppml_cv_pred_dens[i,] <- apply(lik_folds_pred_dens,2,function(x) -log(nrow(lik_folds_pred_dens)) + log(sum(exp(x))))
    apply(lik_folds_pred_dens,2, mean)
    ppml_cv_auc[i,] <- apply(lik_folds_auc,2,mean)
  }
  names(model_fits_cv) <- c('no_offset', 'range_map', 'elevation', 'both_offsets')
  lik_cv_fits <- list(ppml_cv_pred_dens, ppml_cv_auc, model_fits_cv)
  names(lik_cv_fits) <- c('ppml_cv_pred_dens', 'ppml_cv_auc', 'model_fits_cv')
  return(lik_cv_fits)
  
  }



