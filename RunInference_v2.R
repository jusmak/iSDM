RunInference_v2 <- function(wd, species, data, thinning, thinning_PA, weights_area, target_n_obs, n_samples, burn_in) {

  thin_mark <- ifelse(thinning, "thin", "not_thin")
  thin_PA_mark <-  ifelse(thinning_PA, "_thin_PA", "_not_thin_PA")
  weights_mark <- ifelse(weights_area, "wa", "not_wa")
  if(file.exists(paste(wd, "/Model_fits_v2/", species, '_fits_', thin_mark, thin_PA_mark, "_quad_n_",
                       target_n_obs, "_", weights_mark, '.RData', sep = ''))) {
    
    load(paste(wd, "/Model_fits_v2/", species, '_fits_', thin_mark, thin_PA_mark, "_quad_n_",
               target_n_obs, "_", weights_mark, '.RData', sep = ''))
  } else {
    
    source(paste(wd, "R_code/Inference_HB_comb_v2.R", sep = '/'))
    model_fits <- list()
    
    #non-spatial model with single (presence-only) likelihood
    fit_hb <- Inference_non_spat_po(wd, data, n_samples, burn_in)
    model_fits[[1]] <- fit_hb
    
    #non-spatial model with single (presence-absence) likelihood
    fit_hb <- Inference_non_spat_pa(wd, data, n_samples, burn_in)
    model_fits[[2]] <- fit_hb

    #spatial model #1 with single (presence-only) likelihood
    fit_hb <- Inference_spat_po(wd, data, n_samples, burn_in)
    model_fits[[3]] <- fit_hb

    #spatial model #2 with single (presence-only) likelihood
    fit_hb <- Inference_restr_spat_po(wd, data, n_samples, burn_in)
    model_fits[[4]] <- fit_hb
    
    #non-spatial model with combined likelihood
    fit_hb <- Inference_non_spat_comb(wd, data, n_samples, burn_in)
    model_fits[[5]] <- fit_hb

    #spatial model #1
    fit_hb <- Inference_spat_comb(wd, data, n_samples, burn_in)
    model_fits[[6]] <- fit_hb

    #spatial model #2
    fit_hb <- Inference_restr_spat_comb(wd, data, n_samples, burn_in)
    model_fits[[7]] <- fit_hb
    
    names(model_fits) <- c('non_spat_po', 'non_spat_pa', 'spat_po', 'restr_spat_po', 'non_spat_comb', 'spat_comb', 'restr_spat_comb')
    
    save(model_fits, file = paste(wd, "/Model_fits_v2/", species, '_fits_', thin_mark, thin_PA_mark, "_quad_n_",
                                  target_n_obs, "_", weights_mark, '.RData', sep = ''))
  }
  return(model_fits)
}