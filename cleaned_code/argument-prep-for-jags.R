#### WORKFLOW
### 1. Define data to be sent to jags function
### 2. Initialize model parameters
### 3. Define parameters to be tracked
### 4. Define other jags settings

library(bayesm) #to define rwishart()

### 1. Define data to be sent to jags function ----------

#The number of latent classes/ values of true cancer state
nlevel_cancer <- 4
nlevel_cancer_bin <- 2

jags_data_init <-list(
  ## cancer model
  nlevel_cancer=nlevel_cancer, nlevel_cancer_bin=nlevel_cancer_bin, npat=npat,
  cancer_data=cancer_data, npat_cancer_known=npat_cancer_known,
  modmat_cancer=modmat_cancer, npred_cancer=npred_cancer,
  #n_mask=n_mask,
  ## psa model
  nobs_psa=nobs_psa, log_psa_data=log_psa_data, psa_patient_index_map=psa_patient_index_map,
  modmat_ranef_psa=modmat_ranef_psa, modmat_fixef_psa=modmat_fixef_psa, 
  npred_ranef_psa=npred_ranef_psa, npred_fixef_psa=npred_fixef_psa, 
  I_npred_ranef_psa=diag(npred_ranef_psa),
  ## pgg model
  pgg_data=pgg_data, npat_pgg=npat_pgg, pgg_patient_index_map=pgg_patient_index_map,
  modmat_pgg=modmat_pgg, npred_pgg=npred_pgg) 

if(with_mri == "N"){
  jags_data <- jags_data_init
}else{
  if(mri_role == "moderator")){
    jags_data_append <- list(ind_mri_data=ind_mri_data)
  }else if(mri_role == "outcome"){
    jags_data_append <- list(pirads_data = pirads_data, npat_pirads = npat_pirads,
                             pirads_patient_index_map=pirads_patient_index_map,
                             npred_pirads = npred_pirads)
  }else(
    jags_data_append <- list(ind_mri_data=ind_mri_data,
                             pirads_data = pirads_data, npat_pirads = npat_pirads,
                             pirads_patient_index_map=pirads_patient_index_map,
                             npred_pirads = npred_pirads)
  )
}

jags_data <- append(jags_data_init, jags_data_append)



