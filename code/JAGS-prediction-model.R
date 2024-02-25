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
  cancer_int1[i] ~ dnorm(cancer_int1_mean[i], 1) 
  cancer_int2[i] ~ dnorm(cancer_int2_mean[i], 1)
  cancer_int3[i] ~ dnorm(cancer_int3_mean[i], 1) 
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

for (index in 1:npred_ranef_psa[1]) {
	for(k in 1:nlevel_cancer_bin[1]) {
		mu_raw_1[index, k] ~ dnorm(mu_raw[index, k], 0.01)
		mu_raw_2[index, k] ~ dnorm(mu_raw[index, k], 0.01)
		mu_raw_3[index, k] ~ dnorm(mu_raw[index, k], 0.01)
	} 
}

Tau_B_raw_1 ~ dwish(I_npred_ranef_psa[,], (npred_ranef_psa[1]+1))
Sigma_B_raw_1[1:npred_ranef_psa[1], 1:npred_ranef_psa[1]] <- inverse(Tau_B_raw_1[1:npred_ranef_psa[1], 1:npred_ranef_psa[1]])
Tau_B_raw_2 ~ dwish(I_npred_ranef_psa[,], (npred_ranef_psa[2]+1))
Sigma_B_raw_2[1:npred_ranef_psa[2], 1:npred_ranef_psa[2]] <- inverse(Tau_B_raw_2[1:npred_ranef_psa[2], 1:npred_ranef_psa[2]])
Tau_B_raw_3 ~ dwish(I_npred_ranef_psa[,], (npred_ranef_psa[3]+1))
Sigma_B_raw_3[1:npred_ranef_psa[3], 1:npred_ranef_psa[3]] <- inverse(Tau_B_raw_3[1:npred_ranef_psa[3], 1:npred_ranef_psa[3]])

for (i in 1:J){
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
  for(index in 1:npred_pgg[1]) {
    pgg_slope1[index,i] ~ dnorm(pgg_coef_mean[index], 1)
    pgg_slope2[index,i] ~ dnorm(pgg_coef_mean[index], 1)
    pgg_slope3[index,i] ~ dnorm(pgg_coef_mean[index], 1)
  }
  for (index in (npred_pgg[1] + 1):(npred_pgg[1] + (nlevel_cancer[1] - 1))) {
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

## continuation ratio model for true cancer state
for(j in 1:npat[1]){ 
  exponent1_1[j] <- cancer_int1[1] + inprod(cancer_slope1[1:npred_cancer[1],1], modmat_cancer_1[j,1:npred_cancer[1])
  exponent2_1[j] <- cancer_int2[1] + inprod(cancer_slope2[1:npred_cancer[1],1], modmat_cancer_1[j,1:npred_cancer[1]])
  exponent3_1[j] <- cancer_int3[1] + inprod(cancer_slope3[1:npred_cancer[1],1], modmat_cancer_1[j,1:npred_cancer[1]])
  p_eta_1[j, 1] <- exp(exponent1_1[j])/(1+exp(exponent1_1[j]))
  p_eta_1[j, 2] <- 1/(1+exp(exponent1_1[j])) * exp(exponent2_1[j])/(1+exp(exponent2_1[j]))
  p_eta_1[j, 3] <- 1/(1+exp(exponent1_1[j])) * 1/(1+exp(exponent2_1[j])) * exp(exponent3_1[j])/(1+exp(exponent3_1[j]))
  p_eta_1[j, 4] <- 1- p_eta_1[j, 1] - p_eta_1[j, 2] - p_eta_1[j, 3]
}


for(i in 1:npat_cancer_known[1]) {
  cancer_data_1[i] ~ dcat(p_eta_1[i, 1:nlevel_cancer[1]])
  eta_1[i] <- cancer_data_1[i]
  eta.bin_1[i] <- step(eta_1[i]-2) + 1
}

for(i in (npat_cancer_known[1]+1):npat[1]) {
  cancer_state_1[(i-npat_cancer_known[1])] ~ dcat(p_eta_1[i, 1:nlevel_cancer[1]])
  eta_1[i] <- cancer_state_1[(i-npat_cancer_known[1])]
  eta.bin_1[i] <- step(eta_1[i]-2) + 1
}  #for those without SURG

for(j in 1:npat[2]){ 
  exponent1_2[j] <- cancer_int1[2] + inprod(cancer_slope1[1:npred_cancer[2],2], modmat_cancer_2[j,1:npred_cancer[2]])
  exponent2_2[j] <- cancer_int2[2] + inprod(cancer_slope2[1:npred_cancer[2],2], modmat_cancer_2[j,1:npred_cancer[2]])
  exponent3_2[j] <- cancer_int3[2] + inprod(cancer_slope3[1:npred_cancer[2],2], modmat_cancer_2[j,1:npred_cancer[2]])
  p_eta_2[j, 1] <- exp(exponent1_2[j])/(1+exp(exponent1_2[j]))
  p_eta_2[j, 2] <- 1/(1+exp(exponent1_2[j])) * exp(exponent2_2[j])/(1+exp(exponent2_2[j]))
  p_eta_2[j, 3] <- 1/(1+exp(exponent1_2[j])) * 1/(1+exp(exponent2_2[j])) * exp(exponent3_2[j])/(1+exp(exponent3_2[j]))
  p_eta_2[j, 4] <- 1- p_eta_2[j, 1] - p_eta_2[j, 2] - p_eta_2[j, 3]
}

for(i in 1:npat_cancer_known[2]) {
  cancer_data_2[i] ~ dcat(p_eta_2[i, 1:nlevel_cancer[2]])
  eta_2[i] <- cancer_data_2[i]
  eta.bin_2[i] <- step(eta_2[i]-2) + 1
}

for(i in (npat_cancer_known[2]+1):npat[2]) {
  cancer_state_2[(i-npat_cancer_known[2])] ~ dcat(p_eta_2[i, 1:nlevel_cancer[2]])
  eta_2[i] <- cancer_state_2[(i-npat_cancer_known[2])]
  eta.bin_2[i] <- step(eta_2[i]-2) + 1
}  #for those without SURG

for(j in 1:npat[3]){ 
  exponent1_3[j] <- cancer_int1[3] + inprod(cancer_slope1[1:npred_cancer[3],3], modmat_cancer_3[j,1:npred_cancer[3]])
  exponent2_3[j] <- cancer_int2[3] + inprod(cancer_slope2[1:npred_cancer[3],3], modmat_cancer_3[j,1:npred_cancer[3]])
  exponent3_3[j] <- cancer_int3[3] + inprod(cancer_slope3[1:npred_cancer[3],3], modmat_cancer_3[j,1:npred_cancer[3]])
  p_eta_3[j, 1] <- exp(exponent1_3[j])/(1+exp(exponent1_3[j]))
  p_eta_3[j, 2] <- 1/(1+exp(exponent1_3[j])) * exp(exponent2_3[j])/(1+exp(exponent2_3[j]))
  p_eta_3[j, 3] <- 1/(1+exp(exponent1_3[j])) * 1/(1+exp(exponent2_3[j])) * exp(exponent3_3[j])/(1+exp(exponent3_3[j]))
  p_eta_3[j, 4] <- 1- p_eta_3[j, 1] - p_eta_3[j, 2] - p_eta_3[j, 3]
}

for(i in 1:npat_cancer_known[3]) {
  cancer_data_3[i] ~ dcat(p_eta_3[i, 1:nlevel_cancer[3]])
  eta_3[i] <- cancer_data_3[i]
  eta.bin_3[i] <- step(eta_3[i]-2) + 1
}

for(i in (npat_cancer_known[3]+1):npat[3]) {
  cancer_state_3[(i-npat_cancer_known[3])] ~ dcat(p_eta_3[i, 1:nlevel_cancer[3]])
  eta_3[i] <- cancer_state_3[(i-npat_cancer_known[3])]
  eta.bin_3[i] <- step(eta_3[i]-2) + 1
}  #for those without SURG


##linear mixed effects model for PSA
for (i in 1:npat[1]) {
	B_raw_1[i, 1:npred_ranef_psa[1]] ~ dmnorm(mu_raw_1[1:npred_ranef_psa[1], eta.bin_1[i]], Tau_B_raw_1[1:npred_ranef_psa[1], 1:npred_ranef_psa[1]])
	for(index in 1:npred_ranef_psa[1]) {
	  ranef_1[i, index] <- B_raw_1[i, index]
	} 
}

for(j in 1:nobs_psa[1]){
	mu_obs_psa_1[j] <- inprod(ranef_1[psa_patient_index_map_1[j], 1:npred_ranef_psa[1]], modmat_ranef_psa_1[j, 1:npred_ranef_psa[1]]) + inprod(fixef_coefficient[1:npred_fixef_psa[1], 1], modmat_fixef_psa_1[j,1:npred_fixef_psa[1]])
	log_psa_data_1[j] ~ dnorm(mu_obs_psa_1[j], tau_res[1]) 
}

for (i in 1:npat[2]) {
	B_raw_2[i, 1:npred_ranef_psa[2]] ~ dmnorm(mu_raw_2[1:npred_ranef_psa[2], eta.bin_2[i]], Tau_B_raw_2[1:npred_ranef_psa[2], 1:npred_ranef_psa[2]])
	for(index in 1:npred_ranef_psa[2]) {
	  ranef_2[i, index] <- B_raw_2[i, index]
	} 
}

for(j in 1:nobs_psa[2]){
	mu_obs_psa_2[j] <- inprod(ranef_2[psa_patient_index_map_2[j], 1:npred_ranef_psa[2]], modmat_ranef_psa_2[j, 1:npred_ranef_psa[2]]) + inprod(fixef_coefficient[1:npred_fixef_psa[2], 2], modmat_fixef_psa_2[j,1:npred_fixef_psa[2]])
	log_psa_data_2[j] ~ dnorm(mu_obs_psa_2[j], tau_res[2]) 
}

for (i in 1:npat[3]) {
	B_raw_3[i, 1:npred_ranef_psa[3]] ~ dmnorm(mu_raw_3[1:npred_ranef_psa[3], eta.bin_3[i]], Tau_B_raw_3[1:npred_ranef_psa[3], 1:npred_ranef_psa[3]])
	for(index in 1:npred_ranef_psa[3]) {
	  ranef_3[i, index] <- B_raw_3[i, index]
	} 
}

for(j in 1:nobs_psa[3]){
	mu_obs_psa_3[j] <- inprod(ranef_3[psa_patient_index_map_3[j], 1:npred_ranef_psa[3]], modmat_ranef_psa_3[j, 1:npred_ranef_psa[3]]) + inprod(fixef_coefficient[1:npred_fixef_psa[3], 3], modmat_fixef_psa_3[j,1:npred_fixef_psa[3]])
	log_psa_data_3[j] ~ dnorm(mu_obs_psa_3[j], tau_res[3]) 
}

### BIOPSY OUTCOMES AND SURGERY RECEIVED

## continuation ratio model for biopsy upgrading

for(j in 1:npat_pgg[1]){ 
  pgg_exp1_1[j] <- pgg_int1[1] + inprod(pgg_slope1[1:npred_pgg[1],1], modmat_pgg_1[j,1:npred_pgg[1]])+
	               pgg_slope1[(npred_pgg[1] + 1),1] * step(eta_1[pgg_patient_index_map_1[j]] - 2) + 
	               pgg_slope1[(npred_pgg[1] + 2),1] * step(eta_1[pgg_patient_index_map_1[j]] - 3) + 
	               pgg_slope1[(npred_pgg[1] + 3),1] * step(eta_1[pgg_patient_index_map_1[j]] - 4)",
if(mri_role == "MRI_ST"){
  "+ pgg_slope1[(npred_pgg[1] + 4),1] * (pgg_pirads_data_m2_1[pgg_patient_index_map_1[j]]) * (step(eta_1[pgg_patient_index_map_1[j]] - 2)-(1-step(eta_1[pgg_patient_index_map_1[j]] - 2))) * step(pgg_pirads_data_1[pgg_patient_index_map_1[j]] - 3)  
    ## logitP(pgg = 1) w/  (pirads - 2)*[1(eta >= 2)-1(eta=1)]*1(pirads >= 3)"
} else if (mri_role == "MRI"){
  "+ pgg_slope1[(npred_pgg[1] + 4),1] * (step(eta_1[pgg_patient_index_map_1[j]] - 2)-(1-step(eta_1[pgg_patient_index_map_1[j]] - 2))) * step(pgg_pirads_data_1[pgg_patient_index_map_1[j]] - 3)  
    ## logitP(pgg = 1) w/  [1(eta >= 2)-1(eta=1)]*1(pirads >= 3)"
},
"
  pgg_exp2_1[j] <- pgg_int2[1] + inprod(pgg_slope2[1:npred_pgg[1],1], modmat_pgg_1[j,1:npred_pgg[1]])+
	               pgg_slope2[(npred_pgg[1] + 1),1] * step(eta_1[pgg_patient_index_map_1[j]] - 2) + 
	               pgg_slope2[(npred_pgg[1] + 2),1] * step(eta_1[pgg_patient_index_map_1[j]] - 3) + 
	               pgg_slope2[(npred_pgg[1] + 3),1] * step(eta_1[pgg_patient_index_map_1[j]] - 4)",
if(mri_role == "MRI_ST"){
  "+ pgg_slope2[(npred_pgg[1] + 4),1] * (pgg_pirads_data_m2_1[pgg_patient_index_map_1[j]]) * (step(eta_1[pgg_patient_index_map_1[j]] - 2)-(1-step(eta_1[pgg_patient_index_map_1[j]] - 2))) * step(pgg_pirads_data_1[pgg_patient_index_map_1[j]] - 3)  
    ## logitP(pgg = 2|pgg >=2) w/  (pirads - 2)*[1(eta >= 2)-1(eta=1)]*1(pirads >= 3)"
} else if (mri_role == "MRI"){
  "+ pgg_slope2[(npred_pgg[1] + 4),1] * (step(eta_1[pgg_patient_index_map_1[j]] - 2)-(1-step(eta_1[pgg_patient_index_map_1[j]] - 2))) * step(pgg_pirads_data_1[pgg_patient_index_map_1[j]] - 3)  
    ## logitP(pgg = 2|pgg >=2) w/  [1(eta >= 2)-1(eta=1)]*1(pirads >= 3)"
},
"
 pgg_exp3_1[j] <- pgg_int3[1] + inprod(pgg_slope3[1:npred_pgg[1],1], modmat_pgg_1[j,1:npred_pgg[1]])+
	               pgg_slope3[(npred_pgg[1] + 1),1] * step(eta_1[pgg_patient_index_map_1[j]] - 2) + 
	               pgg_slope3[(npred_pgg[1] + 2),1] * step(eta_1[pgg_patient_index_map_1[j]] - 3) + 
	               pgg_slope3[(npred_pgg[1] + 3),1] * step(eta_1[pgg_patient_index_map_1[j]] - 4)",
if(mri_role == "MRI_ST"){
  "+ pgg_slope3[(npred_pgg[1] + 4),1] * (pgg_pirads_data_m2_1[pgg_patient_index_map_1[j]]) * (step(eta_1[pgg_patient_index_map_1[j]] - 2)-(1-step(eta_1[pgg_patient_index_map_1[j]] - 2))) * step(pgg_pirads_data_1[pgg_patient_index_map_1[j]] - 3)  
    ## logitP(pgg = 3|pgg >=3) w/  (pirads - 2)*[1(eta >= 2)-1(eta=1)]*1(pirads >= 3)"
} else if (mri_role == "MRI"){
  "+ pgg_slope3[(npred_pgg[1] + 4),1] * (step(eta_1[pgg_patient_index_map_1[j]] - 2)-(1-step(eta_1[pgg_patient_index_map_1[j]] - 2))) * step(pgg_pirads_data_1[pgg_patient_index_map_1[j]] - 3)  
    ## logitP(pgg = 3|pgg >=3) w/  [1(eta >= 2)-1(eta=1)]*1(pirads >= 3)"
},
"
  p_pgg_1[j, 1] <- exp(pgg_exp1_1[j])/(1+exp(pgg_exp1_1[j]))
  p_pgg_1[j, 2] <- 1/(1+exp(pgg_exp1_1[j])) * exp(pgg_exp2_1[j])/(1+exp(pgg_exp2_1[j]))
  p_pgg_1[j, 3] <- 1/(1+exp(pgg_exp1_1[j])) * 1/(1+exp(pgg_exp2_1[j])) * exp(pgg_exp3_1[j])/(1+exp(pgg_exp3_1[j]))
  p_pgg_1[j, 4] <- 1- p_pgg_1[j, 1] - p_pgg_1[j, 2] - p_pgg_1[j, 3]
}

for(i in 1:npat_pgg[1]) {
  pgg_data_1[i] ~ dcat(p_pgg_1[i,1:nlevel_cancer[1]])
}


for(j in 1:npat_pgg[2]){ 
  pgg_exp1_2[j] <- pgg_int1[2] + inprod(pgg_slope1[1:npred_pgg[1],2], modmat_pgg_2[j,1:npred_pgg[2]])+
	               pgg_slope1[(npred_pgg[2] + 1),2] * step(eta_2[pgg_patient_index_map_2[j]] - 2) + 
	               pgg_slope1[(npred_pgg[2] + 2),2] * step(eta_2[pgg_patient_index_map_2[j]] - 3) + 
	               pgg_slope1[(npred_pgg[2] + 3),2] * step(eta_2[pgg_patient_index_map_2[j]] - 4)",
if(mri_role == "MRI_ST"){
  "+ pgg_slope1[(npred_pgg[2] + 4), 2] * (pgg_pirads_data_m2_2[pgg_patient_index_map_2[j]]) * (step(eta_2[pgg_patient_index_map_2[j]] - 2)-(1-step(eta_2[pgg_patient_index_map_2[j]] - 2))) * step(pgg_pirads_data_2[pgg_patient_index_map_2[j]] - 3)  
    ## logitP(pgg = 1) w/  (pirads - 2)*[1(eta >= 2)-1(eta=1)]*1(pirads >= 3)"
} else if (mri_role == "MRI"){
  "+ pgg_slope1[(npred_pgg[2] + 4), 2] * (step(eta_2[pgg_patient_index_map_2[j]] - 2)-(1-step(eta_2[pgg_patient_index_map_2[j]] - 2))) * step(pgg_pirads_data_2[pgg_patient_index_map_2[j]] - 3)  
    ## logitP(pgg = 1) w/  [1(eta >= 2)-1(eta=1)]*1(pirads >= 3)"
},
"
  pgg_exp2_2[j] <- pgg_int2[2] + inprod(pgg_slope2[1:npred_pgg[2],2], modmat_pgg_2[j,1:npred_pgg[2]])+
	               pgg_slope2[(npred_pgg[2] + 1),2] * step(eta_2[pgg_patient_index_map_2[j]] - 2) + 
	               pgg_slope2[(npred_pgg[2] + 2),2] * step(eta_2[pgg_patient_index_map_2[j]] - 3) + 
	               pgg_slope2[(npred_pgg[2] + 3),2] * step(eta_2[pgg_patient_index_map_2[j]] - 4)",
if(mri_role == "MRI_ST"){
  "+ pgg_slope2[(npred_pgg[2] + 4), 2] * (pgg_pirads_data_m2_2[pgg_patient_index_map_2[j]]) * (step(eta_2[pgg_patient_index_map_2[j]] - 2)-(1-step(eta_2[pgg_patient_index_map_2[j]] - 2))) * step(pgg_pirads_data_2[pgg_patient_index_map_2[j]] - 3)  
    ## logitP(pgg = 2|pgg >=2) w/  (pirads - 2)*[1(eta >= 2)-1(eta=1)]*1(pirads >= 3)"
} else if (mri_role == "MRI"){
  "+ pgg_slope2[(npred_pgg[2] + 4), 2] * (step(eta_2[pgg_patient_index_map_2[j]] - 2)-(1-step(eta_2[pgg_patient_index_map_2[j]] - 2))) * step(pgg_pirads_data_2[pgg_patient_index_map_2[j]] - 3)  
    ## logitP(pgg = 2|pgg >=2) w/  [1(eta >= 2)-1(eta=1)]*1(pirads >= 3)"
},
"
 pgg_exp3_2[j] <- pgg_int3[2] + inprod(pgg_slope3[1:npred_pgg[2],2], modmat_pgg_2[j,1:npred_pgg[2]])+
	               pgg_slope3[(npred_pgg[2] + 1),2] * step(eta_2[pgg_patient_index_map_2[j]] - 2) + 
	               pgg_slope3[(npred_pgg[2] + 2),2] * step(eta_2[pgg_patient_index_map_2[j]] - 3) + 
	               pgg_slope3[(npred_pgg[2] + 3),2] * step(eta_2[pgg_patient_index_map_2[j]] - 4)",
if(mri_role == "MRI_ST"){
  "+ pgg_slope3[(npred_pgg[2] + 4), 2] * (pgg_pirads_data_m2_2[pgg_patient_index_map_2[j]]) * (step(eta_2[pgg_patient_index_map_2[j]] - 2)-(1-step(eta_2[pgg_patient_index_map_2[j]] - 2))) * step(pgg_pirads_data_2[pgg_patient_index_map_2[j]] - 3)  
    ## logitP(pgg =3|pgg>=3) w/  (pirads - 2)*[1(eta >= 2)-1(eta=1)]*1(pirads >= 3)"
} else if (mri_role == "MRI"){
  "+ pgg_slope3[(npred_pgg[2] + 4), 2] * (step(eta_2[pgg_patient_index_map_2[j]] - 2)-(1-step(eta_2[pgg_patient_index_map_2[j]] - 2))) * step(pgg_pirads_data_2[pgg_patient_index_map_2[j]] - 3)  
    ## logitP(pgg =3|pgg>=3) w/  [1(eta >= 2)-1(eta=1)]*1(pirads >= 3)"
},
"
  p_pgg_2[j, 1] <- exp(pgg_exp1_2[j])/(1+exp(pgg_exp1_2[j]))
  p_pgg_2[j, 2] <- 1/(1+exp(pgg_exp1_2[j])) * exp(pgg_exp2_2[j])/(1+exp(pgg_exp2_2[j]))
  p_pgg_2[j, 3] <- 1/(1+exp(pgg_exp1_2[j])) * 1/(1+exp(pgg_exp2_2[j])) * exp(pgg_exp3_2[j])/(1+exp(pgg_exp3_2[j]))
  p_pgg_2[j, 4] <- 1- p_pgg_2[j, 1] - p_pgg_2[j, 2] - p_pgg_2[j, 3]
}

for(i in 1:npat_pgg[2]) {
  pgg_data_2[i] ~ dcat(p_pgg_2[i,1:nlevel_cancer[2]])
}

for(j in 1:npat_pgg[3]){ 
  pgg_exp1_3[j] <- pgg_int1[3] + inprod(pgg_slope1[1:npred_pgg[3],3], modmat_pgg_3[j,1:npred_pgg[3]])+
	               pgg_slope1[(npred_pgg[3] + 1),3] * step(eta_3[pgg_patient_index_map_3[j]] - 2) + 
	               pgg_slope1[(npred_pgg[3] + 2),3] * step(eta_3[pgg_patient_index_map_3[j]] - 3) + 
	               pgg_slope1[(npred_pgg[3] + 3),3] * step(eta_3[pgg_patient_index_map_3[j]] - 4)",
if(mri_role == "MRI_ST"){
  "+ pgg_slope1[(npred_pgg[3] + 4),3] * (pgg_pirads_data_m2_3[pgg_patient_index_map_3[j]]) * (step(eta_3[pgg_patient_index_map_3[j]] - 2)-(1-step(eta_3[pgg_patient_index_map_3[j]] - 2))) * step(pgg_pirads_data_3[pgg_patient_index_map_3[j]] - 3)  
    ## logitP(pgg = 1) w/  (pirads - 2)*[1(eta >= 2)-1(eta=1)]*1(pirads >= 3)"
} else if (mri_role == "MRI"){
  "+ pgg_slope1[(npred_pgg[3] + 4),3] * (step(eta_3[pgg_patient_index_map_3[j]] - 2)-(1-step(eta_3[pgg_patient_index_map_3[j]] - 2))) * step(pgg_pirads_data_3[pgg_patient_index_map_3[j]] - 3)  
    ## logitP(pgg = 1) w/  [1(eta >= 2)-1(eta=1)]*1(pirads >= 3)"
},
"
  pgg_exp2_1[j] <- pgg_int2[3] + inprod(pgg_slope2[1:npred_pgg[3],3], modmat_pgg_3[j,1:npred_pgg[3]])+
	               pgg_slope2[(npred_pgg[3] + 1),3] * step(eta_3[pgg_patient_index_map_3[j]] - 2) + 
	               pgg_slope2[(npred_pgg[3] + 2),3] * step(eta_3[pgg_patient_index_map_3[j]] - 3) + 
	               pgg_slope2[(npred_pgg[3] + 3),3] * step(eta_3[pgg_patient_index_map_3[j]] - 4)",
if(mri_role == "MRI_ST"){
  "+ pgg_slope2[(npred_pgg[3] + 4),3] * (pgg_pirads_data_m2_3[pgg_patient_index_map_3[j]]) * (step(eta_3[pgg_patient_index_map_3[j]] - 2)-(1-step(eta_3[pgg_patient_index_map_3[j]] - 2))) * step(pgg_pirads_data_3[pgg_patient_index_map_3[j]] - 3)  
    ## logitP(pgg = 2|pgg >=2) w/ (pirads - 2)*[1(eta >= 3)-1(eta<=2)]*1(pirads >= 3)"
} else if (mri_role == "MRI"){
  "+ pgg_slope2[(npred_pgg[3] + 4),3] * (step(eta_3[pgg_patient_index_map_3[j]] - 2)-(1-step(eta_3[pgg_patient_index_map_3[j]] - 2))) * step(pgg_pirads_data_3[pgg_patient_index_map_3[j]] - 3)  
    ## logitP(pgg = 2|pgg >=2) w/ [1(eta >= 3)-1(eta<=2)]*1(pirads >= 3)"
},
"
 pgg_exp3_3[j] <- pgg_int3[3] + inprod(pgg_slope3[1:npred_pgg[3],3], modmat_pgg_3[j,1:npred_pgg[3]])+
	               pgg_slope3[(npred_pgg[3] + 1),3] * step(eta_3[pgg_patient_index_map_3[j]] - 2) + 
	               pgg_slope3[(npred_pgg[3] + 2),3] * step(eta_3[pgg_patient_index_map_3[j]] - 3) + 
	               pgg_slope3[(npred_pgg[3] + 3),3] * step(eta_3[pgg_patient_index_map_3[j]] - 4)",
if(mri_role == "MRI_ST"){
  "+ pgg_slope3[(npred_pgg[3] + 4),3] * (pgg_pirads_data_m2_3[pgg_patient_index_map_3[j]]) * (step(eta_3[pgg_patient_index_map_3[j]] - 2)-(1-step(eta_3[pgg_patient_index_map_3[j]] - 2))) * step(pgg_pirads_data_3[pgg_patient_index_map_3[j]] - 3)  
    ## logitP(pgg =3|pgg>=3) w/ (pirads - 2)*[1(eta >= 3)-1(eta<=2)]*1(pirads >= 3)"
} else if (mri_role == "MRI"){
  "+ pgg_slope3[(npred_pgg[3] + 4),3] * (step(eta_3[pgg_patient_index_map_3[j]] - 2)-(1-step(eta_3[pgg_patient_index_map_3[j]] - 2))) * step(pgg_pirads_data_3[pgg_patient_index_map_3[j]] - 3)  
    ## logitP(pgg =3|pgg>=3) w/ [1(eta >= 3)-1(eta<=2)]*1(pirads >= 3)"
},
"
  p_pgg_3[j, 1] <- exp(pgg_exp1_3[j])/(1+exp(pgg_exp1_3[j]))
  p_pgg_3[j, 2] <- 1/(1+exp(pgg_exp1_3[j])) * exp(pgg_exp2_3[j])/(1+exp(pgg_exp2_3[j]))
  p_pgg_3[j, 3] <- 1/(1+exp(pgg_exp1_3[j])) * 1/(1+exp(pgg_exp2_3[j])) * exp(pgg_exp3_3[j])/(1+exp(pgg_exp3_3[j]))
  p_pgg_3[j, 4] <- 1- p_pgg_3[j, 1] - p_pgg_3[j, 2] - p_pgg_3[j, 3]
}

for(i in 1:npat_pgg[3]) {
  pgg_data_3[i] ~ dcat(p_pgg_3[i,1:nlevel_cancer[3]])
}

",

"}",
  file = paste(location.of.r.scripts, 
                jags_file_name, sep="/"))


