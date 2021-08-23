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
K <- 4
K.bin <- 2

jags_data<-list(K=K, K.bin=K.bin, n=n,
                eta.data=eta.data, n_eta_known=n_eta_known,
                V.ETA=V.ETA.data, d.V.ETA=d.V.ETA,

                n_obs_psa=n_obs_psa, Y=Y, subj_psa=subj_psa,
                Z=Z.data, X=X.data, d.Z=d.Z, d.X=d.X, I_d.Z=diag(d.Z),

               # NPC=NPC, n_npc=n_npc, subj_npc=subj_npc,
                #V.NPC=V.NPC.data, d.V.NPC=d.V.NPC, #NCS.offset=NCS.offset,

                PGG=PGG, n_pgg=n_pgg, subj_pgg=subj_pgg,
                V.PGG=V.PGG.data, d.V.PGG=d.V.PGG#,

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

rho_int0 <- rnorm((K-1),0,1)
rho_coef <- rnorm(d.V.ETA, mean=0, sd=0.25)

eta.hat <- pt.data$bx.pgg[is.na(eta.data)]

xi <- rlnorm(d.Z)
mu_raw <- as.matrix(cbind(rnorm(d.Z),rnorm(d.Z)))
Tau_B_raw <- rwishart((d.Z+1),diag(d.Z)*var_vec)$W
sigma_res <- min(rlnorm(1),1)

beta <- rnorm(d.X)

#nu.BX <- rnorm((d.U.BX+1), mean=0, sd=0.1)

#beta.NPC <- rnorm((d.V.NPC + 1), mean=0, sd=0.5 )
#beta.NPC[(d.V.NPC+1)] <- abs(beta.NPC[(d.V.NPC+1)])
#r.NPC <- runif(1,1,10)

alpha0 <- rnorm((K-1),0,1)
gamma.PGG <- rnorm((d.V.PGG+(K-1)),mean=0,sd=0.25)
gamma.PGG[(d.V.PGG+1):(d.V.PGG+(K-1))] <- abs(gamma.PGG[(d.V.PGG+1):(d.V.PGG+(K-1))])

#beta.MPC <- rnorm((d.V.MPC + 1), mean=0, sd=0.5 )
#beta.MPC[(d.V.MPC+1)] <- abs(beta.MPC[(d.V.MPC+1)])
#r.MPC <- runif(1,1,10)

#beta.LAT <- rnorm((d.V.LAT + 1), mean=0, sd=0.25 )
#beta.LAT[(d.V.LAT+1)] <- abs(beta.LAT[(d.V.LAT+1)])

#omega.SURG <-c (rnorm((d.W.SURG+2), mean=0, sd=0.01))

list(rho_int0=rho_int0, rho_coef=rho_coef, #cat_int=cat_int,
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
params <- c("rho_int", "rho_coef",
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
