Set_data_spde = function(wd, species, thinning, target_n_obs, weights_area, data) {
  
  #read data if a data file has been created earlier
  thin_mark = ifelse(thinning, 'thin', 'not_thin')
  weights_mark = ifelse(weights_area, 'wa', 'not_wa')
  
  if (file.exists(paste0('Compiled_data/Inference_list_', species, '_', thin_mark, '_quad_n_',
                       target_n_obs, '_', weights_mark, '.R'))) {
    
    load(paste0('Compiled_data/Inference_list_', species, '_', thin_mark, '_quad_n_',
               target_n_obs, '_', weights_mark, '.R'))
    
  } else {  
  
    #for INLA, create a mesh for approximating the Gaussian random field
    #create a mesh
    hull = inla.nonconvex.hull(data$PO_training_data$coordinates, convex = -0.01, resolution = 500)
    #### check how to define a mesh size  ####
    mesh = inla.mesh.2d(boundary = hull, cutoff = 20, offset = c(40, 80), max.edge = c(160, 320))
    
    #cross-validation folds - only complete checklists
    n_folds = 4
    
    #interpolation
    set.seed(1)
    fold_po = sample(1:n_folds,size=sum(data$PA_training_data$response==1),replace=TRUE)
    fold_ab = sample(1:n_folds,size=sum(data$PA_training_data$response==0),replace=TRUE)
    
    list_ind = function(ind, fold_po, fold_ab, response_data) {
      ind = sort(c(which(response_data==1)[fold_po == ind], which(response_data==0)[fold_ab == ind]))
      return(ind)
    }
    
    int_folds = lapply(1:n_folds, function(x) list_ind(x, fold_po, fold_ab, data$PA_training_data$response))
    
    
    #extrapolation - 20 stripes and set to four folds
    fold_po = kmeans(data$PA_training_data$coordinates, n_folds*5)
    fold_po = fold_po$cluster
    
    n_runs = 20000
    temp_runs = matrix(NA, n_runs, 1)
    for (i in 1:n_runs) {
      set.seed(i)
      ind_runs = sample(rep(1:n_folds, 5),20,replace = F)
      temp = lapply(1:n_folds, function(x) which(ind_runs==x))
      prev_temp = lapply(temp, function(x) which(fold_po %in% x))
      prev_temp = lapply(prev_temp, function(x) sum(data$PA_training_data$response[x]))
      temp_runs[i] = sd(unlist(prev_temp))
    }
    
    set.seed(which.min(temp_runs))
    ind_runs = sample(rep(1:n_folds, 5),20,replace = F)
    temp = lapply(1:n_folds, function(x) which(ind_runs==x))
    ext_folds = lapply(temp, function(x) which(fold_po %in% x))
    
    names(int_folds) = paste0('Interpolation_set', 1:4)
    names(ext_folds) = paste0('Extrapolation_set', 1:4)
    pa_testing_folds = c(int_folds, ext_folds)
    
    
    #set training and testing set indices for presence-only data. This is for the second extrapolation test.
    #find presence points that are located in the testing area of the presence-absence data
    #set extent of the testing presence-absence points
    ind_pres = which(data$PO_training_data$response == 1)
    temp_backgr = (max(ind_pres)+1):length(data$PO_training_data$response)

    po_training_folds = vector(mode = 'list', length = n_folds)
    
    for (k in 1:length(po_training_folds)) {
      d_mat = matrix(NA, nrow = length(data$PA_training_data$response[ext_folds[[k]]]), ncol = length(ind_pres))
      for (i in 1:nrow(d_mat)) {
        for (j in 1:ncol(d_mat)) {
          d_mat[i,j] = sqrt((data$PA_training_data$coordinates[ext_folds[[k]][i],1]-data$PO_training_data$coordinates[ind_pres[j],1])^2+
                               (data$PA_training_data$coordinates[ext_folds[[k]][i],2]-data$PO_training_data$coordinates[ind_pres[j],2])^2)
        }
      }
      #find cells which are more than 30 km away from the presence-absence points
      ind_keep = which(apply(d_mat, 2, function(x) any(x<30))==F)
      po_training_folds[[k]] = c(ind_pres[ind_keep], temp_backgr)
    }
    
    po_training_folds = append(rep(list(1:length(data$PO_training_data$response)),4), po_training_folds)
    names(po_training_folds) = names(pa_testing_folds)
    
    
    data_inla = list('po_coordinates' = data$PO_training_data$coordinates, 'po_covariates' = data$PO_training_data$covariates,
                      'po_response' = data$PO_training_data$response, 'po_weights' = data$PO_training_data$weights,
                      'po_offset' = data$PO_training_data$offset_expert,
                      'pa_coordinates' = data$PA_training_data$coordinates,
                      'pa_covariates' = data$PA_training_data$covariates, 'pa_response' = data$PA_training_data$response,
                      'pa_offset' = data$PA_training_data$offset_expert,
                      'cov_mean' = data$PO_training_data$cov_mean, 'cov_sd' = data$PO_training_data$cov_sd,
                      'spde_hull' = hull, 'spde_mesh' = mesh, 'pa_testing_folds' = pa_testing_folds,
                      'po_training_folds' = po_training_folds)    
    

    save(data_inla, file = paste0('Compiled_data/Inference_list_', species, '_', thin_mark, '_quad_n_',
                                 target_n_obs, '_', weights_mark, '.R'))    
  }
  
  return(data_inla)

}