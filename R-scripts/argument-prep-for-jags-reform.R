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

jags_data<-list(K=nlevel_cancer, K.bin=nlevel_cancer_bin, n=npat,
                eta.data=cancer_data, n_eta_known=npat_cancer_known,
                V.ETA=modmat_cancer, d.V.ETA=npred_cancer,

                n_obs_psa=nobs_psa, Y=log_psa_data, subj_psa=psa_patient_index_map,
                Z=modmat_ranef_psa, X=modmat_fixef_psa, d.Z=npred_ranef_psa, d.X=npred_fixef_psa, 
                I_d.Z=diag(npred_ranef_psa),

               # NPC=NPC, n_npc=n_npc, subj_npc=subj_npc,
                #V.NPC=V.NPC.data, d.V.NPC=d.V.NPC, #NCS.offset=NCS.offset,

                PGG=pgg_data, n_pgg=npat_pgg, subj_pgg=pgg_patient_index_map,
                V.PGG=modmat_pgg, d.V.PGG=npred_pgg#,

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

eta.hat <- pt.data$bx.pgg[is.na(cancer_data)]

xi <- rlnorm(npred_ranef_psa)
mu_raw <- as.matrix(cbind(rnorm(npred_ranef_psa),rnorm(npred_ranef_psa)))
Tau_B_raw <- rwishart((npred_ranef_psa+1),diag(npred_ranef_psa)*var_vec)$W
sigma_res <- min(rlnorm(1),1)

beta <- rnorm(npred_fixef_psa)

#nu.BX <- rnorm((d.U.BX+1), mean=0, sd=0.1)

#beta.NPC <- rnorm((d.V.NPC + 1), mean=0, sd=0.5 )
#beta.NPC[(d.V.NPC+1)] <- abs(beta.NPC[(d.V.NPC+1)])
#r.NPC <- runif(1,1,10)

alpha0 <- rnorm((nlevel_cancer-1),0,1)
gamma.PGG <- rnorm((npred_pgg+(nlevel_cancer-1)),mean=0,sd=0.25)
gamma.PGG[(npred_pgg+1):(npred_pgg+(nlevel_cancer-1))] <- abs(gamma.PGG[(npred_pgg+1):(npred_pgg+(nlevel_cancer-1))])

#beta.MPC <- rnorm((d.V.MPC + 1), mean=0, sd=0.5 )
#beta.MPC[(d.V.MPC+1)] <- abs(beta.MPC[(d.V.MPC+1)])
#r.MPC <- runif(1,1,10)

#beta.LAT <- rnorm((d.V.LAT + 1), mean=0, sd=0.25 )
#beta.LAT[(d.V.LAT+1)] <- abs(beta.LAT[(d.V.LAT+1)])

#omega.SURG <-c (rnorm((d.W.SURG+2), mean=0, sd=0.01))

list(cancer_intercept0=cancer_intercept0, cancer_slope=cancer_slope, #cat_int=cat_int,
     eta.hat=eta.hat,
     xi=xi, mu_raw=mu_raw, Tau_B_raw=Tau_B_raw, sigma_res=sigma_res, beta=beta,
 #    beta.NPC=beta.NPC, r.NPC=r.NPC,
     alpha0=alpha0, gamma.PGG=gamma.PGG#,
#     beta.MPC=beta.MPC, r.MPC=r.MPC,
#     beta.LAT=beta.LAT#,
 #    nu.BX=nu.BX, omega.SURG=omega.SURG
     ) }
#




### 3. Define parameters to be tracked
params <- c("cancer_intercept", "cancer_slope",
            "eta.hat",
            "xi", "mu_int", "mu_slope",
            "sigma_int", "sigma_slope", "sigma_res",  "cov_int_slope",
            "b.vec",
            "beta",
#            "beta.NPC", "r.NPC", #"mu_npc",
            "alpha", "gamma.PGG"#, #"p_rc",
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






