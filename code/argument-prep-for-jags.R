########### Argumennt prep for a hierarchical model using J sets of synthetic data. 


#### WORKFLOW
### 1. Define data to be sent to jags function
### 2. Initialize model parameters
### 3. Define parameters to be tracked
### 4. Define other jags settings

library(bayesm) #to define rwishart()

### 1. Define data to be sent to jags function ----------

# The number of latent classes/ values of true cancer state
nlevel_cancer <- rep(4,3)
nlevel_cancer_bin <- rep(2,3)
J <- 3 #number of institutions for hierarchical model

jags_data_init <-list(
  J = J,
  ## prior means
  cancer_int1_mean = cancer_int1_mean, cancer_int2_mean = cancer_int2_mean, cancer_int3_mean = cancer_int3_mean,
  pgg_int1_mean = pgg_int1_mean, pgg_int2_mean = pgg_int2_mean, pgg_int3_mean = pgg_int3_mean,
  
  ## cancer model
  nlevel_cancer=nlevel_cancer, nlevel_cancer_bin=nlevel_cancer_bin, npat=npat,
  cancer_data_1=cancer_data_true[[1]], npat_cancer_known=npat_cancer_known,
  modmat_cancer_1 = modmat_cancer[[1]], npred_cancer=npred_cancer,
  cancer_data_2=cancer_data_true[[2]], modmat_cancer_2 = modmat_cancer[[2]],
  cancer_data_3=cancer_data_true[[3]], modmat_cancer_3 = modmat_cancer[[3]],
  
  
  ## psa model
  nobs_psa=nobs_psa, log_psa_data_1=log_psa_data[[1]], psa_patient_index_map_1=psa_patient_index_map[[1]],
  modmat_ranef_psa_1=modmat_ranef_psa[[1]], modmat_fixef_psa_1=modmat_fixef_psa[[1]], 
  npred_ranef_psa=npred_ranef_psa, npred_fixef_psa=npred_fixef_psa, 
  I_npred_ranef_psa=diag(npred_ranef_psa[1]),
  log_psa_data_2=log_psa_data[[2]], psa_patient_index_map_2=psa_patient_index_map[[2]],
  modmat_ranef_psa_2=modmat_ranef_psa[[2]], modmat_fixef_psa_2=modmat_fixef_psa[[2]], 
  log_psa_data_3=log_psa_data[[3]], psa_patient_index_map_3=psa_patient_index_map[[3]],
  modmat_ranef_psa_3=modmat_ranef_psa[[3]], modmat_fixef_psa_3=modmat_fixef_psa[[3]], 
  
  ## biopsy model
  pgg_data_1=pgg_data[[1]], npat_pgg=npat_pgg, pgg_patient_index_map_1=pgg_patient_index_map[[1]],
  modmat_pgg_1=modmat_pgg[[1]], npred_pgg=npred_pgg, 
  pgg_data_2=pgg_data[[2]], pgg_patient_index_map_2=pgg_patient_index_map[[2]],
  modmat_pgg_2=modmat_pgg[[2]], 
  pgg_data_3=pgg_data[[3]], pgg_patient_index_map_3=pgg_patient_index_map[[3]],
  modmat_pgg_3=modmat_pgg[[3]]) 

if(mri_role == 0){
  jags_data <- jags_data_init
}else{
  if(mri_role == "MRI"){
    jags_data_append <- list(pgg_pirads_data_1 = pgg_pirads_data[[1]],
                             pgg_pirads_data_2 = pgg_pirads_data[[2]],
                             pgg_pirads_data_3 = pgg_pirads_data[[3]],
                             pgg_pirads_data_m2_1 = pgg_pirads_data_m2[[1]],
                             pgg_pirads_data_m2_2 = pgg_pirads_data_m2[[2]],
                             pgg_pirads_data_m2_3 = pgg_pirads_data_m2[[3]])
  }
  jags_data <- append(jags_data_init, jags_data_append)
}


### 2. Initialize model parameters ---------------
inits <- function(){
  cancer_state_1 <- cancer_state[[1]]$pgg_mean
  cancer_state_2 <- cancer_state[[2]]$pgg_mean
  cancer_state_3 <- cancer_state[[3]]$pgg_mean
  #latent cancer model
  cancer_int1 <- cancer_int2 <- cancer_int3 <- rnorm(J,0,1)
  cancer_slope1 <- matrix(rnorm(npred_cancer[1]*J,mean=0,sd=0.25), ncol = J)
  cancer_slope2 <- matrix(rnorm(npred_cancer[1]*J,mean=0,sd=0.25), ncol = J)
  cancer_slope3 <- matrix(rnorm(npred_cancer[1]*J,mean=0,sd=0.25), ncol = J)
  cancer_coef_mean = rnorm(npred_cancer[1],mean=0,sd=0.25)
  
  #psa model
  mu_raw_1 <- as.matrix(cbind(rnorm(npred_ranef_psa[1]),rnorm(npred_ranef_psa[1])))
  Tau_B_raw_1 <- rwishart((npred_ranef_psa[1]+1),diag(npred_ranef_psa[1])*var_vec[[1]])$W
  resid_var_psa <- rep(min(rlnorm(1),1), J)
  fixef_coefficient <- matrix(rnorm(npred_fixef_psa[1] * J), ncol = J)
  
  mu_raw_2 <- as.matrix(cbind(rnorm(npred_ranef_psa[2]),rnorm(npred_ranef_psa[2])))
  Tau_B_raw_2 <- rwishart((npred_ranef_psa[2]+1),diag(npred_ranef_psa[2])*var_vec[[2]])$W

  mu_raw_3 <- as.matrix(cbind(rnorm(npred_ranef_psa[3]),rnorm(npred_ranef_psa[3])))
  Tau_B_raw_3 <- rwishart((npred_ranef_psa[3]+1),diag(npred_ranef_psa[3])*var_vec[[3]])$W

  # Biopsy model
  pgg_int1 <- pgg_int2 <-pgg_int3 <-rnorm(J,0,1)
  if(mri_role==0){
    pgg_coef_mean = rnorm(npred_pgg[1] + (nlevel_cancer[1]-1),mean=0,sd=0.25)
    pgg_slope1 <- matrix(rnorm(J*(npred_pgg[1] + (nlevel_cancer[1]-1)),mean=0,sd=0.25), ncol = J)
    pgg_slope2 <- matrix(rnorm(J*(npred_pgg[1] + (nlevel_cancer[1]-1)),mean=0,sd=0.25), ncol = J)
    pgg_slope3 <- matrix(rnorm(J*(npred_pgg[1] + (nlevel_cancer[1]-1)),mean=0,sd=0.25), ncol = J)
  }else{
    pgg_coef_mean = rnorm(npred_pgg[1] + (nlevel_cancer[1]),mean=0,sd=0.25)
    pgg_slope1 <- matrix(rnorm(J*(npred_pgg[1] + nlevel_cancer[1]),mean=0,sd=0.25), ncol = J)
    pgg_slope2 <- matrix(rnorm(J*(npred_pgg[1] + nlevel_cancer[1]),mean=0,sd=0.25), ncol = J)
    pgg_slope3 <- matrix(rnorm(J*(npred_pgg[1] + nlevel_cancer[1]),mean=0,sd=0.25), ncol = J)
  }

  
  inits_list <- list(cancer_state_1 = cancer_state_1, cancer_state_2 = cancer_state_2, cancer_state_3 = cancer_state_3,
                     cancer_int1=cancer_int1, cancer_int2=cancer_int2, cancer_int3=cancer_int3, 
                     cancer_slope1=cancer_slope1,cancer_slope2=cancer_slope2,cancer_slope3=cancer_slope3,
                     cancer_coef_mean = cancer_coef_mean,
                     
                     mu_raw_1=mu_raw_1, mu_raw_2=mu_raw_2, mu_raw_3=mu_raw_3, 
                     Tau_B_raw_1=Tau_B_raw_1, Tau_B_raw_2=Tau_B_raw_2, Tau_B_raw_3=Tau_B_raw_3, 
                     resid_var_psa=resid_var_psa, 
                     fixef_coefficient=fixef_coefficient, 
                     
                     pgg_int1=pgg_int1, pgg_int2=pgg_int2, pgg_int3=pgg_int3, 
                     pgg_slope1=pgg_slope1,pgg_slope2=pgg_slope2,pgg_slope3=pgg_slope3,
                     pgg_coef_mean = pgg_coef_mean)
  
  inits_list
}

### 3. Define parameters to be tracked ---------------
params <- c(
            "cancer_int1", "cancer_int2", "cancer_int3",
            "cancer_slope1", "cancer_slope2", "cancer_slope3", 
            "eta_1", "eta_2", "eta_3",
            "p_eta_1", "p_eta_2", "p_eta_3",
            
            "pgg_int1", "pgg_int2", "pgg_int3",
            "pgg_slope1", "pgg_slope2", "pgg_slope3"
)



### 4. Define other jags settings -------------------- 
#### it is possible for chain to converge with many fewer iterations. Check the appropriate number of iterations for your dataset.
#### change length; burn-in; number thinned; number of chains
## n.iter 10k minimum
## n.burnin 10k minimum
## n.chains == # cores - 1

# change length; burn-in; number thinned; number of chains
n.iter <- 10000; n.burnin <- 2500; n.thin <- 10; n.chains <- 1






