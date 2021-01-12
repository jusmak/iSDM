Inference_ppml_comb <- function(wd, training_data, PA_data, offset_full) {
  
  
  form = ~ cov1 + cov2 + cov3 + I(cov1^2) + I(cov2^2) + I(cov3^2) + 
    I(cov4^2) + I(cov1*cov2) + I(cov1*cov3) + I(cov1*cov4) + I(cov2*cov3) +
    I(cov2*cov4) + I(cov3*cov4)
  
  bias_form_PO = ~ cov4
  bias_form_PA = ~ cov4
  

  #set data as a data frame in the ppmlasso format
  #coordinates, covariates, response and weights
  #keep PO data first in the order
  
  ppm_PO_data <- data.frame(cbind(training_data$coordinates, 
                                  training_data$covariates,
                   training_data$response,training_data$weights))
  
  ppm_PA_data <- data.frame(cbind(PA_data$coordinates, PA_data$covariates,
                                  PA_data$response,rep(1,nrow(PA_data$covariates))))
  
  colnames(ppm_PO_data) <- c('X', 'Y', colnames(training_data$covariates),
                             'Pres', 'wt')
  colnames(ppm_PA_data) <- c('X', 'Y', colnames(training_data$covariates),
                             'Pres', 'wt')
  
  #fit_ppml = ppmlasso(form, data = ppm_PO_data, n.fits = 50, max.it = 10)
  
  source(paste(wd, "/R_code/Renner_et_al_2019_comb_lasso_copy_JM.R", sep = ''))
  
  #check link function
  fit_comb_ppml <- comb_lasso(env_formula = form,
                              bias_formula = list(bias_form_PO, bias_form_PA),
                              intercept_env = 1,
                              intercept_bias = 1,
                              quad_data = ppm_PO_data[ppm_PO_data$Pres==0,1:6],
                              sp_data = list(ppm_PO_data[ppm_PO_data$Pres==1,1:6],
                                             ppm_PA_data[,1:6]),
                              sp_y = list(rep(1,sum(ppm_PO_data$Pres)), ppm_PA_data$Pres),
                              dat.type = c("PO", "OCC"),
                              coord = c('X','Y'),
                              sp_res = 1,
                              link = "logit",
                              site.area = 1,
                              standardise = TRUE,
                              ob_temp = ppm_PO_data$wt)
  
  #store model formula, coefficients, likelihood and BIC
  beta_temp <- fit_comb_ppml$beta[1:(length(fit_comb_ppml$beta)-2)]
  names(fit_comb_ppml$beta) <- NULL
  fit_comb_temp <- list(fit_comb_ppml$formula[[1]], fit_comb_ppml$beta, fit_comb_ppml$pen_likelihood,
                        fit_comb_ppml$criterion.matrix[,3], fit_comb_ppml$penalty_vec)
  names(fit_comb_temp) <- c('formula', 'coefficients', 'pen_lik', 'BIC', 'pen_vector')
  return(fit_comb_temp)
}