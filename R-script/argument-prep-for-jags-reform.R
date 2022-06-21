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



### 1. Define data to be sent to jags function

#The number of latent classes/ values of true cancer state
nlevel_cancer <- 4
nlevel_cancer_bin <- 2

# jags_data<-list(K=K, K.bin=K.bin, n=n,
#                 eta.data=eta.data, n_eta_known=n_eta_known,
#                 V.ETA=V.ETA.data, d.V.ETA=d.V.ETA,
#                 
#                 n_obs_psa=n_obs_psa, Y=Y, subj_psa=subj_psa,
#                 Z=Z.data, X=X.data, d.Z=d.Z, d.X=d.X, I_d.Z=diag(d.Z),
# 
#                 PGG=PGG, n_pgg=n_pgg, subj_pgg=subj_pgg,
#                 V.PGG=V.PGG.data, d.V.PGG=d.V.PGG#,
# ) 

jags_data<-list(nlevel_cancer=nlevel_cancer, nlevel_cancer_bin=nlevel_cancer_bin, npat=npat,
                cancer_data=cancer_data, npat_cancer_known=npat_cancer_known,
                modmat_cancer=modmat_cancer, npred_cancer=npred_cancer,

                nobs_psa=nobs_psa, log_psa_data=log_psa_data, psa_patient_index_map=psa_patient_index_map,
                modmat_ranef_psa=modmat_ranef_psa, modmat_fixef_psa=modmat_fixef_psa, 
                npred_ranef_psa=npred_ranef_psa, npred_fixef_psa=npred_fixef_psa, 
                I_npred_ranef_psa=diag(npred_ranef_psa),

               # NPC=NPC, n_npc=n_npc, subj_npc=subj_npc,
                #V.NPC=V.NPC.data, d.V.NPC=d.V.NPC, #NCS.offset=NCS.offset,

               pgg_data=pgg_data, npat_pgg=npat_pgg, pgg_patient_index_map=pgg_patient_index_map,
               modmat_pgg=modmat_pgg, npred_pgg=npred_pgg#,

                #MPC=MPC, n_mpc=n_mpc, subj_mpc=subj_mpc,
                #V.MPC=V.MPC.data, d.V.MPC=d.V.MPC,

                #LAT=LAT, n_lat=n_lat, subj_lat=subj_lat,
                #V.LAT=V.LAT.data, d.V.LAT=d.V.LAT#,

                #n_bx=n_bx, BX=BX, subj_bx=subj_bx,
                #U.BX=U.BX.data, d.U.BX=d.U.BX ,
                #n_surg=n_surg, SURG=SURG, subj_surg=subj_surg,
                #W.SURG=W.SURG.data, d.W.SURG=d.W.SURG

                ) #,



### 2. Initialize model parameters

inits <- function() {

cancer_intercept0 <- rnorm((nlevel_cancer-1),0,1)
cancer_slope <- rnorm(npred_cancer, mean=0, sd=0.25)

cancer_state <- pt.data$bx.pgg[is.na(cancer_data)]

scale_ranef_mean_psa <- rlnorm(npred_ranef_psa)
mu_raw <- as.matrix(cbind(rnorm(npred_ranef_psa),rnorm(npred_ranef_psa)))
Tau_B_raw <- rwishart((npred_ranef_psa+1),diag(npred_ranef_psa)*var_vec)$W
resid_var_psa <- min(rlnorm(1),1)

fixef_coefficient <- rnorm(npred_fixef_psa)


pgg_intercept0 <- rnorm((nlevel_cancer-1),0,1)
pgg_slope <- rnorm((npred_pgg+(nlevel_cancer-1)),mean=0,sd=0.25)
pgg_slope[(npred_pgg+1):(npred_pgg+(nlevel_cancer-1))] <- abs(pgg_slope[(npred_pgg+1):(npred_pgg+(nlevel_cancer-1))])


list(cancer_intercept0=cancer_intercept0, cancer_slope=cancer_slope, #cat_int=cat_int,
     cancer_state=cancer_state,
     scale_ranef_mean_psa=scale_ranef_mean_psa, mu_raw=mu_raw, Tau_B_raw=Tau_B_raw, 
     resid_var_psa=resid_var_psa, fixef_coefficient=fixef_coefficient,
 #    beta.NPC=beta.NPC, r.NPC=r.NPC,
 pgg_intercept0=pgg_intercept0, pgg_slope=pgg_slope#,
#     beta.MPC=beta.MPC, r.MPC=r.MPC,
#     beta.LAT=beta.LAT#,
 #    nu.BX=nu.BX, omega.SURG=omega.SURG
     ) }
#




### 3. Define parameters to be tracked
#reformat Yate's code
## rho_int -> cancer_intercept
## rho_coef -> cancer_slope
## eta.hat -> cancer_state
## xi -> scale_ranef_mean_psa
## mu_int -> ranef_intercept; mu_slope  -> ranef_slope (mean of)
## sigma_int -> ranef_var_intercept; sigma_slope -> ranef_var_slope; cov_int_slope -> ranef_cov
## sigma_res ->  resid_var_psa
## b.vec -> ranef
## beta -> fixef_coefficient
## alpha -> pgg_intercept
## gamma.PGG -> pgg_slope
params <- c("cancer_intercept", "cancer_slope",
            "cancer_state",
            "scale_ranef_mean_psa", "ranef_intercept", "ranef_slope",
            "ranef_var_intercept", "ranef_var_slope", "ranef_cov","resid_var_psa",  
            "ranef",
            "fixef_coefficient",
#            "beta.NPC", "r.NPC", #"mu_npc",
            "pgg_intercept", "pgg_slope"#, #"p_rc",
#            "beta.MPC", "r.MPC", #"mu_mpc",
#            "beta.LAT"#, #"p_lat"
            #"nu.BX", #"p_bx",
            #"omega.SURG"#, #"p_surg"
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






