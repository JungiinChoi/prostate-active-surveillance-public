#Yates Coley
#rycoley@gmail.com
#Last UPDATED July 13, 2017
#Annotations updated December 21, 2017
#This script defines arguments needed to run JAGS


#IOP components are commented out because they do not improve model estimation
#NPC, MPC, and LAT components are commented out because they do not improve model estimation.




#### WORKFLOW
### 1. Define data to be sent to jags function
### 2. Initialize model parameters
### 3. Define parameters to be tracked
### 4. Define other jags settings

library(bayesm) #to define rwishart()

### 1. Define data to be sent to jags function

#The number of latent classes/ values of true cancer state
nlevel_cancer <- 4
nlevel_cancer_bin <- 2

jags_data<-list(nlevel_cancer=nlevel_cancer, nlevel_cancer_bin=nlevel_cancer_bin, npat=npat,
                cancer_data=cancer_data, npat_cancer_known=npat_cancer_known,
                modmat_cancer=modmat_cancer, npred_cancer=npred_cancer,
                alpha = 1, beta = 1, 
                n_mask=n_mask,
                nobs_psa=nobs_psa, log_psa_data=log_psa_data, psa_patient_index_map=psa_patient_index_map,
                modmat_ranef_psa=modmat_ranef_psa, modmat_fixef_psa=modmat_fixef_psa, 
                npred_ranef_psa=npred_ranef_psa, npred_fixef_psa=npred_fixef_psa, 
                I_npred_ranef_psa=diag(npred_ranef_psa),

                # # NPC=NPC, n_npc=n_npc, subj_npc=subj_npc,
                # #V.NPC=V.NPC.data, d.V.NPC=d.V.NPC, #NCS.offset=NCS.offset,
                
                # pgg_data=pgg_data, npat_pgg=npat_pgg, pgg_patient_index_map=pgg_patient_index_map,
                # modmat_pgg=modmat_pgg, npred_pgg=npred_pgg,
                
                pirads_data = pirads_data, npat_pirads = npat_pirads,pirads_patient_index_map=pirads_patient_index_map,
                npred_pirads = npred_pirads,
                
                abserror_data=abserror_data, npat_abserror=npat_abserror, 
                modmat_abserror=modmat_abserror, npred_abserror=npred_abserror
               #  
               #  layer_eta2_data = layer_eta2_data, layer_eta3_data=layer_eta3_data,
               # # npat_layer = npat_layer, 
               #  npat_layer_eta2 = npat_layer_eta2, npat_layer_eta3 = npat_layer_eta3,
               # # layer_patient_index_map = layer_patient_index_map,
               #  # layer_eta2_patient_index_map = layer_eta2_patient_index_map,
               #  # layer_eta3_patient_index_map = layer_eta3_patient_index_map,
               #  layer_state_index = layer_state_index,
               #  layer_state_map = layer_state_map
              
                
) #,




### 2. Initialize model parameters

inits <- function() {
  
  cancer_int1 <- cancer_int2 <-cancer_int3 <-rnorm(1,0,1)
  cancer_slope1 <- rnorm(npred_cancer,mean=0,sd=0.25)
  cancer_slope2 <- rnorm(npred_cancer,mean=0,sd=0.25)
  cancer_slope3 <- rnorm(npred_cancer,mean=0,sd=0.25)
  cancer_coef_mean = rnorm(npred_cancer,mean=0,sd=0.25)
  cancer_inv_sigma_sq = rgamma(1, 1, 1)
  
  cancer_state <- pt.data$bx.pgg[is.na(cancer_data)]         ## set initial values for cancer state to be biopsy value     
  abserror_state <- pt.data$diff_eta_pgg[is.na(abserror_data)] ## initial values all NA, choose from prior
  layer_state <- rep(NA_real_, length(layer_state_index))
  
  scale_ranef_mean_psa <- rlnorm(npred_ranef_psa)
  mu_raw <- as.matrix(cbind(rnorm(npred_ranef_psa),rnorm(npred_ranef_psa)))
  Tau_B_raw <- rwishart((npred_ranef_psa+1),diag(npred_ranef_psa)*var_vec)$W
  resid_var_psa <- min(rlnorm(1),1)
  
  fixef_coefficient <- rnorm(npred_fixef_psa)
  
  
  # pgg_int1 <- pgg_int2 <-pgg_int3 <-rnorm(1,0,1)
  # pgg_slope1 <- rnorm(npred_pgg + (nlevel_cancer-1),mean=0,sd=0.25)
  # pgg_slope2 <- rnorm(npred_pgg + (nlevel_cancer-1),mean=0,sd=0.25)
  # pgg_slope3 <- rnorm(npred_pgg + (nlevel_cancer-1),mean=0,sd=0.25)
  # pgg_coef_mean = rnorm(npred_pgg + (nlevel_cancer-1),mean=0,sd=0.25)
  # pgg_coef_mean[(npred_pgg+1):(npred_pgg+(nlevel_cancer-1))] = abs(pgg_coef_mean[(npred_pgg+1):(npred_pgg+(nlevel_cancer-1))])
  # pgg_inv_sigma_sq = rgamma(1, 1, 1)
  # pgg_slope1[(npred_pgg+1):(npred_pgg+(nlevel_cancer-1))] <- abs(pgg_slope1[(npred_pgg+1):(npred_pgg+(nlevel_cancer-1))])
  # pgg_slope2[(npred_pgg+1):(npred_pgg+(nlevel_cancer-1))] <- abs(pgg_slope2[(npred_pgg+1):(npred_pgg+(nlevel_cancer-1))])
  # pgg_slope3[(npred_pgg+1):(npred_pgg+(nlevel_cancer-1))] <- abs(pgg_slope3[(npred_pgg+1):(npred_pgg+(nlevel_cancer-1))])
  # # 
  pirads_int1 <- pirads_int2 <- abs(rnorm(1,0,1))
  pirads_inv_sigma_sq = rgamma(1, 1, 1)
  pirads_slope1 <- rnorm(npred_pirads + (nlevel_cancer-1),mean=0,sd=0.25)
  pirads_slope2 <- rnorm(npred_pirads + (nlevel_cancer-1),mean=0,sd=0.25)
  pirads_coef_mean = rnorm(npred_pirads + (nlevel_cancer-1),mean=0,sd=0.25) 
  pirads_coef_mean[(npred_pirads+1):(npred_pirads+(nlevel_cancer-1))] = abs(pirads_coef_mean[(npred_pirads+1):(npred_pirads+(nlevel_cancer-1))])
  pirads_slope1[(npred_pirads+1):(npred_pirads+(nlevel_cancer-1))] = abs(pirads_slope1[(npred_pirads+1):(npred_pirads+(nlevel_cancer-1))])
  pirads_slope2[(npred_pirads+1):(npred_pirads+(nlevel_cancer-1))] = abs(pirads_slope2[(npred_pirads+1):(npred_pirads+(nlevel_cancer-1))])
  #
  abserror_int1 <- abserror_int2 <-abserror_int3 <-rnorm(1,0,1)
  abserror_slope <- rnorm(npred_abserror + (nlevel_cancer-1),mean=0,sd=0.25)
  abserror_coef_mean = rnorm(npred_abserror+ (nlevel_cancer-1),mean=0,sd=0.25)
  abserror_inv_sigma_sq = rgamma(1, 1, 1)
  # 
  #layer_eta2_int = layer_eta3_int = abs(rnorm(1,0,1))
  #
  
  list(cancer_int1=cancer_int1, cancer_int2=cancer_int2, cancer_int3=cancer_int3, 
       cancer_slope1=cancer_slope1,cancer_slope2=cancer_slope2,cancer_slope3=cancer_slope3,
       cancer_inv_sigma_sq = cancer_inv_sigma_sq, 
       cancer_coef_mean = cancer_coef_mean,
       
       scale_ranef_mean_psa=scale_ranef_mean_psa, mu_raw=mu_raw, Tau_B_raw=Tau_B_raw, 
       resid_var_psa=resid_var_psa, fixef_coefficient=fixef_coefficient,
       # 
       # pgg_int1=pgg_int1, pgg_int2=pgg_int2, pgg_int3=pgg_int3, 
       # pgg_slope1=pgg_slope1,pgg_slope2=pgg_slope2,pgg_slope3=pgg_slope3,
       # pgg_inv_sigma_sq = pgg_inv_sigma_sq, 
       # pgg_coef_mean = pgg_coef_mean,
       
       pirads_int1 = pirads_int1, pirads_int2= pirads_int2,
       pirads_slope1=pirads_slope1,pirads_slope2=pirads_slope2,
       pirads_inv_sigma_sq = pirads_inv_sigma_sq, 
       pirads_coef_mean = pirads_coef_mean,
       
       abserror_int1 = abserror_int1, abserror_int2 = abserror_int2, abserror_int3 = abserror_int3,
       abserror_slope = abserror_slope,
       abserror_coef_mean = abserror_coef_mean,
       abserror_inv_sigma_sq = abserror_inv_sigma_sq
       
       # layer_eta2_int = layer_eta2_int, layer_eta3_int = layer_eta3_int,
       # 
       # cancer_state=cancer_state,
       # abserror_state = abserror_state,
       # layer_state=layer_state
       
  ) }
# 




### 3. Define parameters to be tracked


params <- c("eta.track",
            "cancer_int1", "cancer_int2", "cancer_int3",
            "cancer_slope1", "cancer_slope2", "cancer_slope3", 
            "cancer_sigma_sq",
            "cancer_state",
            
            "scale_ranef_mean_psa", "ranef_intercept", "ranef_slope",
            "ranef_var_intercept", "ranef_var_slope", "ranef_cov","resid_var_psa",  
            "ranef",
            "fixef_coefficient",
            
            # "pgg_int1", "pgg_int2", "pgg_int3",
            # "pgg_slope1", "pgg_slope2", "pgg_slope3", 
            # "pgg_sigma_sq",
            
            "pirads_int1", "pirads_int2", 
            "pirads_slope1", "pirads_slope2", "pirads_slope3", 
            "pirads_sigma_sq",
            
            "abserror_int1", "abserror_int2", "abserror_int3",
            "abserror_slope",
            "abserror_coef_mean",
            "abserror_inv_sigma_sq",
            "abserror_state"
            
            # "layer_eta2_int" ,
            # "layer_eta3_int",
            # "layer_state" 
)



### 4. Define other jags settings
#### it is possible for chain to converge with many fewer iterations. Check the appropriate number of iterations for your dataset.
#### change length; burn-in; number thinned; number of chains
## n.iter 10k minimum
## n.burnin 10k minimum
## n.chains == # cores - 1

# change length; burn-in; number thinned; number of chains
n.iter <- 10000; n.burnin <- 2500; n.thin <- 10; n.chains <- 1
#n.iter <- 40000; n.burnin <- 10000; n.thin <- 20; n.chains <- 1
#n.iter <- 100; n.burnin <- 20; n.thin <- 5; n.chains <- 1

#### I run 1 chain here and run 5 models in parallel (instead of 5 chains back-to-back)


