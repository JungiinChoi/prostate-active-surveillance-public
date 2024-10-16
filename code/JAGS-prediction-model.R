if(K>1){
  jags_file_name <- paste0("cv-JAGS-prediction-model-",mri_role,seed,".txt")
}else{
  jags_file_name <- paste0("JAGS-prediction-model-",mri_role, seed,".txt")
}


cat(
  "model {
  
### PRIORS FOR LATENT CANCER STATE MODEL

# priors for sequential model for ordinal eta as outcome
for (index in 1:npred_cancer[1]) {
   cancer_coef_mean[index] ~ dnorm(0,1)
}
  
for (i in 1:J){
  for(j in 1:max(npat)){
    cancer_int1[i,j] ~ dnorm(cancer_int1_mean[i], 30) 
    cancer_int2[i,j] ~ dnorm(cancer_int2_mean[i], 30)
    cancer_int3[i,j] ~ dnorm(cancer_int3_mean[i], 30) 
  }
  for(index in 1:npred_cancer[i]) {
    cancer_slope1[index,i] ~ dnorm(cancer_coef_mean[index], 1)
    cancer_slope2[index,i] ~ dnorm(cancer_coef_mean[index], 1)
    cancer_slope3[index,i] ~ dnorm(cancer_coef_mean[index], 1)
  }
}


### PRIORS FOR PSA MIXED MODEL

for (index in 1:npred_ranef_psa[1]) {
	for(k in 1:nlevel_cancer_bin[1]) {
		mu_raw[index, k] ~ dnorm(0, 1)
	} 
}

for (i in 1:J){
  for (index in 1:npred_ranef_psa[i]) {
	  for(k in 1:nlevel_cancer_bin[i]) {
		  mu_raw_i[i, index, k] ~ dnorm(mu_raw[index, k], 0.01)
	  }
  }
  Tau_B_raw[i,1:2,1:2] ~ dwish(I_npred_ranef_psa[,], (npred_ranef_psa[i]+1))
  resid_var_psa[i] ~ dunif(0, 1)
  tau_res[i] <- pow(resid_var_psa[i], -2)
  
  for(index in 1:npred_fixef_psa[i]) {
  	fixef_coefficient[index, i] ~ dnorm(0, 0.01)
  }
}

### PRIORS FOR BIOPSY MODEL

for (index in 1:npred_pgg[1]) {
  pgg_coef_mean[index] ~ dnorm(0,1)
}
for (index in (npred_pgg[1] + 1):(npred_pgg[1] + (nlevel_cancer[1] - 1))) {
  pgg_coef_mean[index] ~ dnorm(0,1)
}
",
if(mri_role %in% c("MRI", "MRI_ST")){"pgg_coef_mean[npred_pgg[1] + nlevel_cancer[1]] ~ dnorm(0,1)"},

"\n
for (i in 1:J){
  pgg_int1[i] ~ dnorm(pgg_int1_mean[i], 1) 
  pgg_int2[i] ~ dnorm(pgg_int2_mean[i], 1)
  pgg_int3[i] ~ dnorm(pgg_int3_mean[i], 1)
  for(index in 1:npred_pgg[i]) {
    pgg_slope1[index,i] ~ dnorm(pgg_coef_mean[index], 1)
    pgg_slope2[index,i] ~ dnorm(pgg_coef_mean[index], 1)
    pgg_slope3[index,i] ~ dnorm(pgg_coef_mean[index], 1)
  }
  for (index in (npred_pgg[i] + 1):(npred_pgg[i] + (nlevel_cancer[i] - 1))) {
    pgg_slope1[index,i] ~ dnorm(pgg_coef_mean[index], 1)
    pgg_slope2[index,i] ~ dnorm(pgg_coef_mean[index], 1)
    pgg_slope3[index,i] ~ dnorm(pgg_coef_mean[index], 1)
  }
}
",

if(mri_role %in% c("MRI", "MRI_ST")){
  "
  for (i in 1:J){
    pgg_slope1[npred_pgg[1] + nlevel_cancer[1],i] ~ dnorm(pgg_coef_mean[npred_pgg[1] + nlevel_cancer[1]],1)
    pgg_slope2[npred_pgg[1] + nlevel_cancer[1],i] ~ dnorm(pgg_coef_mean[npred_pgg[1] + nlevel_cancer[1]],1)
    pgg_slope3[npred_pgg[1] + nlevel_cancer[1],i] ~ dnorm(pgg_coef_mean[npred_pgg[1] + nlevel_cancer[1]],1)
  }"
},

"

### LIKELIHOOD


for (i in 1:J){
  ## continuation ratio model for true cancer state
  for(j in 1:npat[i]){ 
    exponent1[i,j] <- cancer_int1[i,j] + inprod(cancer_slope1[1:npred_cancer[i],i], modmat_cancer[i,j,1:npred_cancer[i]])
    exponent2[i,j] <- cancer_int2[i,j] + inprod(cancer_slope2[1:npred_cancer[i],i], modmat_cancer[i,j,1:npred_cancer[i]])
    exponent3[i,j] <- cancer_int3[i,j] + inprod(cancer_slope3[1:npred_cancer[i],i], modmat_cancer[i,j,1:npred_cancer[i]])
    p_eta[i, j, 1] <- exp(exponent1[i,j])/(1+exp(exponent1[i,j]))
    p_eta[i, j, 2] <- 1/(1+exp(exponent1[i,j])) * exp(exponent2[i,j])/(1+exp(exponent2[i,j]))
    p_eta[i, j, 3] <- 1/(1+exp(exponent1[i,j])) * 1/(1+exp(exponent2[i,j])) * exp(exponent3[i,j])/(1+exp(exponent3[i,j]))
    p_eta[i, j, 4] <- 1- p_eta[i, j, 1] - p_eta[i, j, 2] - p_eta[i, j, 3]
  }

  for (j in 1:npat_cancer_known[i]){
    cancer_data[i,j] ~ dcat(p_eta[i, j, 1:nlevel_cancer[i]])
    eta[i,j] <- cancer_data[i,j]
    eta.bin[i,j] <- step(eta[i,j]-2) + 1
  }
  
  for(j in (npat_cancer_known[i]+1):npat[i]) {
    cancer_state[i,j-npat_cancer_known[i]] ~ dcat(p_eta[i,j, 1:4])
    eta[i,j] <- cancer_state[i,(j-npat_cancer_known[i])]
    eta.bin[i,j] <- step(eta[i,j]-2) + 1
  }  #for those without SURG
  for (j in (npat[i]+1):max(npat)){
    eta[i,j] <- 0
  }
  for (j in (npat[i]-npat_cancer_known[i]+1):max(npat - npat_cancer_known)){
    cancer_state[i,j] ~ dcat(c(1/4,1/4,1/4,1/4))
  }
}


##linear mixed effects model for PSA

for (i in 1:J){
  for (j in 1:npat[i]) {
  	B_raw[i, j, 1:npred_ranef_psa[i]] ~ dmnorm(mu_raw_i[i,1:npred_ranef_psa[1], eta.bin[i,j]], Tau_B_raw[i,1:npred_ranef_psa[i], 1:npred_ranef_psa[i]])
  	for(index in 1:npred_ranef_psa[i]) {
  	  ranef[i, j, index] <- B_raw[i,j, index]
  	} 
  }
  for(j in 1:nobs_psa[i]){
  	mu_obs_psa[i,j] <- inprod(ranef[i, psa_patient_index_map[i,j], 1:npred_ranef_psa[i]], modmat_ranef_psa[i, j, 1:npred_ranef_psa[i]]) + inprod(fixef_coefficient[1:npred_fixef_psa[i], i], modmat_fixef_psa[i,j,1:npred_fixef_psa[i]])
  	log_psa_data[i,j] ~ dnorm(mu_obs_psa[i,j], tau_res[i]) 
  }
}

### BIOPSY OUTCOMES AND SURGERY RECEIVED

## continuation ratio model for biopsy upgrading

for (i in 1:J){
  for(j in 1:npat_pgg[i]){ 
    pgg_exp1[i,j] <- pgg_int1[i] + inprod(pgg_slope1[1:npred_pgg[i],1], modmat_pgg[i,j,1:npred_pgg[i]])+
  	               pgg_slope1[(npred_pgg[i] + 1),1] * step(eta[i, pgg_patient_index_map[i,j]] - 2) + 
  	               pgg_slope1[(npred_pgg[i] + 2),1] * step(eta[i, pgg_patient_index_map[i,j]] - 3) + 
  	               pgg_slope1[(npred_pgg[i] + 3),1] * step(eta[i, pgg_patient_index_map[i,j]] - 4)",
  
  if(mri_role == "MRI_ST"){
    "+ pgg_slope1[(npred_pgg[i] + 4),i] * (pgg_pirads_data_m2[i,pgg_patient_index_map[i,j]]) * (step(eta[i, pgg_patient_index_map[i,j]] - 2)-(1-step(eta[1, pgg_patient_index_map[i,j]] - 2))) * step(pgg_pirads_data[i,pgg_patient_index_map[i,j]] - 3)  
      ## logitP(pgg = 1) w/  (pirads - 2)*[1(eta >= 2)-1(eta=1)]*1(pirads >= 3)"
  } else if (mri_role == "MRI"){
    "+ pgg_slope1[(npred_pgg[i] + 4),i] * (step(eta[i, pgg_patient_index_map[i,j]] - 2)-(1-step(eta[1, pgg_patient_index_map[i,j]] - 2))) * step(pgg_pirads_data[i,pgg_patient_index_map[i,j]] - 3)  
      ## logitP(pgg = 1) w/  [1(eta >= 2)-1(eta=1)]*1(pirads >= 3)"
  },

  "
    pgg_exp2[i,j] <- pgg_int2[i] + inprod(pgg_slope2[1:npred_pgg[i],2], modmat_pgg[i,j,1:npred_pgg[i]])+
  	               pgg_slope2[(npred_pgg[i] + 1),2] * step(eta[i, pgg_patient_index_map[i,j]] - 2) + 
  	               pgg_slope2[(npred_pgg[i] + 2),2] * step(eta[i, pgg_patient_index_map[i,j]] - 3) + 
  	               pgg_slope2[(npred_pgg[i] + 3),2] * step(eta[i, pgg_patient_index_map[i,j]] - 4)",
  
  if(mri_role == "MRI_ST"){
    "+ pgg_slope2[(npred_pgg[i] + 4),i] * (pgg_pirads_data_m2[i,pgg_patient_index_map[i,j]]) * (step(eta[i, pgg_patient_index_map[i,j]] - 2)-(1-step(eta[1, pgg_patient_index_map[i,j]] - 2))) * step(pgg_pirads_data[i,pgg_patient_index_map[i,j]] - 3)  
      ## logitP(pgg = 1) w/  (pirads - 2)*[1(eta >= 2)-1(eta=1)]*1(pirads >= 3)"
  } else if (mri_role == "MRI"){
    "+ pgg_slope2[(npred_pgg[i] + 4),i] * (step(eta[i, pgg_patient_index_map[i,j]] - 2)-(1-step(eta[1, pgg_patient_index_map[i,j]] - 2))) * step(pgg_pirads_data[i,pgg_patient_index_map[i,j]] - 3)  
      ## logitP(pgg = 1) w/  [1(eta >= 2)-1(eta=1)]*1(pirads >= 3)"
  },

  "
    pgg_exp3[i,j] <- pgg_int3[i] + inprod(pgg_slope3[1:npred_pgg[i],3], modmat_pgg[i,j,1:npred_pgg[i]])+
  	               pgg_slope3[(npred_pgg[i] + 1),3] * step(eta[i, pgg_patient_index_map[i,j]] - 2) + 
  	               pgg_slope3[(npred_pgg[i] + 2),3] * step(eta[i, pgg_patient_index_map[i,j]] - 3) + 
  	               pgg_slope3[(npred_pgg[i] + 3),3] * step(eta[i, pgg_patient_index_map[i,j]] - 4)",
  
  if(mri_role == "MRI_ST"){
    "+ pgg_slope3[(npred_pgg[i] + 4),i] * (pgg_pirads_data_m2[i,pgg_patient_index_map[i,j]]) * (step(eta[i, pgg_patient_index_map[i,j]] - 2)-(1-step(eta[1, pgg_patient_index_map[i,j]] - 2))) * step(pgg_pirads_data[i,pgg_patient_index_map[i,j]] - 3)  
      ## logitP(pgg = 1) w/  (pirads - 2)*[1(eta >= 2)-1(eta=1)]*1(pirads >= 3)"
  } else if (mri_role == "MRI"){
    "+ pgg_slope3[(npred_pgg[i] + 4),i] * (step(eta[i, pgg_patient_index_map[i,j]] - 2)-(1-step(eta[1, pgg_patient_index_map[i,j]] - 2))) * step(pgg_pirads_data[i,pgg_patient_index_map[i,j]] - 3)  
      ## logitP(pgg = 1) w/  [1(eta >= 2)-1(eta=1)]*1(pirads >= 3)"
  },

  "
    logit(p_pgg[i, j, 1]) <- max(0.0001, pgg_exp1[i,j])
    p_pgg[i, j, 2] <- max(0.0001, 1/(1+exp(pgg_exp1[i,j])) * exp(pgg_exp2[i,j])/(1+exp(pgg_exp2[i,j])))
    p_pgg[i, j, 3] <- max(0.0001, 1/(1+exp(pgg_exp1[i,j])) * 1/(1+exp(pgg_exp2[i,j])) * exp(pgg_exp3[i,j])/(1+exp(pgg_exp3[i,j])))
    p_pgg[i, j, 4] <- max(0.0001, 1- p_pgg[i, j, 1] - p_pgg[i, j, 2] - p_pgg[i, j, 3])
  }
  
  for(j in 1:npat_pgg[i]) {
    pgg_data[i,j] ~ dcat(p_pgg[i,j,1:nlevel_cancer[i]])
  }

}

",

"}",
  file = paste(location.of.r.scripts, 
                jags_file_name, sep="/"))


