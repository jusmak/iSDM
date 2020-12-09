Inference_glmnet_cv_x_fold <- function(wd, training_data, offset_full, ind_fold) {
  
  glmnet_cv_pred_dens <- matrix(NA, ncol = length(offset_full), nrow = length(offset_full))
  glmnet_cv_auc <- matrix(NA, ncol = length(offset_full), nrow = length(offset_full))
  model_fits_cv <- list()
  
  for (i in 1:length(offset_full)) {
    model_fits_cv_offset <- list()
    lik_folds_pred_dens <- matrix(NA, nrow = length(ind_fold), ncol = length(offset_full))
    lik_folds_auc <- matrix(NA, nrow = length(ind_fold), ncol = length(offset_full))
    
    for (j in 1:length(ind_fold)) {
      
      ind_temp <- ind_fold[[j]]
      ind_temp <- (1:sum(training_data$response==1))[-ind_temp]
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
      offset_temp <- offset_full[[i]][ind_keep_1,]
      pop_temp <- matrix(apply(offset_temp, 2, function(x) sum(x*weights_temp)), nrow = 1)
      offset_temp <- mapply(function(x,y) x*sum(resp_temp)/y, data.frame(offset_temp), data.frame(pop_temp))
      
      offset_temp_1 = log(offset_temp[,1]) - mean(log(offset_temp[,1])) + 
        log(offset_temp[,2]) - mean(log(offset_temp[,2]))

      form = ~ poly(cov1, degree = 2, raw = TRUE) + poly(cov2, degree = 2, raw = TRUE) + 
        poly(cov3, degree = 2, raw = TRUE) + poly(cov4, degree = 2, raw = TRUE)
      cv.glmnet.1 <-glmnet(formula = form, x = cov_temp_sc, 
                              y = resp_temp/weights_temp,
                              family = "poisson",
                              offset = offset_temp_1,
                              weights = weights_temp)

      model_fits_cv_offset[[j]] <- cv.glmnet.1

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
        
        #get offset of the closest training point
        offset_temp <- offset_full[[k]][ind_keep_1,]
        pop_temp <- matrix(apply(offset_temp, 2, function(x) sum(x*weights_temp)), nrow = 1)
        offset_temp <- mapply(function(x,y) x*sum(resp_temp)/y, data.frame(offset_temp), data.frame(pop_temp))
        
        #offset for the quadrature points
        offset_temp_test <- offset_temp[resp_temp==0,]
        
        #offset for the presence points
        ind_close <- matrix(NA, nrow = sum(resp_temp_test==1))
        for (m in 1:length(ind_close)) {
          ind_close[m] <- which.min((coord_temp_test[m,1] - coord_temp[,1])^2 +
                                      (coord_temp_test[m,2] - coord_temp[,2])^2)
        }
        
        offset_temp_test <- rbind(offset_temp[ind_close,], offset_temp_test)
        
        pred_lik <- matrix(NA, nrow = nrow(cov_temp_test_sc), ncol = 1)
        offset_temp <- log(offset_temp_test[,1]) - mean(log(offset_temp_test[,1])) + 
          log(offset_temp_test[,2]) - mean(log(offset_temp_test[,2]))
        
        pred_validation <- predict(object = cv.glmnet.1, newx = cov_temp_test_sc,
                newoffset = offset_temp)

        #predictive density
        pred_lik <- dpois(resp_temp_test, exp(pred_validation[,ncol(pred_validation)])*weights_temp_test, log = TRUE)
        lik_folds_pred_dens[j,k] <- sum(pred_lik)
        
        #AUC
        occ_prob <- 1-exp(-exp(pred_validation[,ncol(pred_validation)] + log(weights_temp_test)))
        suppressMessages(
          myRoc <- pROC::roc(response = resp_temp_test, predictor = occ_prob))
        lik_folds_auc[j,k] <- myRoc$auc
      }
    }
    model_fits_cv[[i]] <- model_fits_cv_offset
    glmnet_cv_pred_dens[i,] <- apply(lik_folds_pred_dens,2,mean)
    glmnet_cv_auc[i,] <- apply(lik_folds_auc,2,mean)
  }
  names(model_fits_cv) <- c('no_offset', 'range_map', 'elevation', 'both_offsets')
  lik_cv_fits <- list(glmnet_cv_pred_dens, glmnet_cv_auc, model_fits_cv)
  names(lik_cv_fits) <- c('glmnet_cv_pred_dens', 'glmnet_cv_auc', 'model_fits_cv')
  return(lik_cv_fits)
  
  }



