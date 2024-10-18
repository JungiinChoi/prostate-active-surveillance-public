########### Argumennt prep for a hierarchical model using J sets of synthetic data. 


#### WORKFLOW
### 1. Define data to be sent to jags function
### 2. Initialize model parameters
### 3. Define parameters to be tracked
### 4. Define other jags settings

library(bayesm) #to define rwishart()

### 1. Define data to be sent to jags function ----------

# The number of latent classes/ values of true cancer state
J <- 3 #number of institutions for hierarchical model
nlevel_cancer <- rep(4,J)
nlevel_cancer_bin <- rep(2,J)

jags_data_init <-list(
  J = J,
  ## prior means
  cancer_int1_mean = cancer_int1_mean, cancer_int2_mean = cancer_int2_mean, cancer_int3_mean = cancer_int3_mean,
  pgg_int1_mean = pgg_int1_mean, pgg_int2_mean = pgg_int2_mean, pgg_int3_mean = pgg_int3_mean,
  
  ## cancer model
  nlevel_cancer=nlevel_cancer, nlevel_cancer_bin=nlevel_cancer_bin, npat=as.integer(npat),
  cancer_data = cancer_data_true, modmat_cancer = modmat_cancer,
  npat_cancer_known=npat_cancer_known, npred_cancer=npred_cancer,
  
  ## psa model
  nobs_psa=nobs_psa, log_psa_data=log_psa_data, psa_patient_index_map=psa_patient_index_map,
  modmat_ranef_psa=modmat_ranef_psa, modmat_fixef_psa=modmat_fixef_psa, 
  npred_ranef_psa=npred_ranef_psa, npred_fixef_psa=npred_fixef_psa, 
  I_npred_ranef_psa=diag(npred_ranef_psa[1]),
  
  ## biopsy model
  pgg_data = pgg_data, npat_pgg=npat_pgg, pgg_patient_index_map=pgg_patient_index_map,
  npred_pgg=npred_pgg, modmat_pgg = modmat_pgg)

if(mri_role == 0){
  jags_data <- jags_data_init
}else{
  if(mri_role == "MRI" | mri_role == "MRI_ST"){
    jags_data_append <- list(pgg_pirads_data = pgg_pirads_data,
                             pgg_pirads_data_m2 = pgg_pirads_data_m2)
  }
  jags_data <- append(jags_data_init, jags_data_append)
}


### 2. Initialize model parameters ---------------
inits <- function(){
  cancer_state <- cancer_state
  #latent cancer model
  cancer_int1 <- cancer_int2 <- cancer_int3 <- matrix(rnorm(max(npat)*J,0,1),nrow = J)
  cancer_slope1 <- matrix(rnorm(npred_cancer[1]*J,mean=0,sd=0.25), ncol = J)
  cancer_slope2 <- matrix(rnorm(npred_cancer[1]*J,mean=0,sd=0.25), ncol = J)
  cancer_slope3 <- matrix(rnorm(npred_cancer[1]*J,mean=0,sd=0.25), ncol = J)
  cancer_coef_mean = rnorm(npred_cancer[1],mean=0,sd=0.25)
  
  #psa model
  mu_raw_i <- array(0, dim = c(J,2,2))
  for (i in 1:J){
    mu_raw_i[i,,] <- cbind(rnorm(npred_ranef_psa[i]),rnorm(npred_ranef_psa[i]))
  }
  
  Tau_B_raw <- array(0, dim = c(J,npred_ranef_psa[1],npred_ranef_psa[1]))
  for (i in 1:J){
    Tau_B_raw[i,,] <- rwishart((npred_ranef_psa[i]+1),diag(npred_ranef_psa[i])*var_vec[[i]])$W
  }
  resid_var_psa <- rep(min(rlnorm(1),1), J)
  fixef_coefficient <- matrix(rnorm(npred_fixef_psa[1] * J), ncol = J)

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

  
  inits_list <- list(cancer_state = cancer_state, 
                     cancer_int1=cancer_int1, cancer_int2=cancer_int2, cancer_int3=cancer_int3, 
                     cancer_slope1=cancer_slope1,cancer_slope2=cancer_slope2,cancer_slope3=cancer_slope3,
                     cancer_coef_mean = cancer_coef_mean,
                     
                     mu_raw_i=mu_raw_i,
                     Tau_B_raw=Tau_B_raw,
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
            "eta", "p_eta",
            
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
n.iter <- 1000; n.burnin <- 250; n.thin <- 1; n.chains <- 1






