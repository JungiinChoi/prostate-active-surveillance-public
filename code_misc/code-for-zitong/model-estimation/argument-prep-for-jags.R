#Yates Coley
#rycoley@gmail.com
#2018 August 23
#This script defines arguments needed to run JAGS


#### WORKFLOW
### 0. Load R packages
### 1. Define data to be sent to jags function
### 2. Initialize model parameters
### 3. Define parameters to be tracked



library(bayesm) #to define rwishart()


### 1. Define data to be sent to jags function

#The number of latent classes/ values of true cancer state
K <- 4 
K.bin <- 2

jags_data<-list(K=K, K.bin=K.bin, n=n, 
                eta.data=eta.data, n_eta_known=n_eta_known, 
                V.ETA=V.ETA.data, d.V.ETA=d.V.ETA,
                
                n_obs_psa=n_obs_psa, Y=Y, subj_psa=subj_psa, 
                Z=Z.data, X=X.data, d.Z=d.Z, d.X=d.X, I_d.Z=diag(d.Z), 
                
                PGG=PGG, n_pgg=n_pgg, subj_pgg=subj_pgg, 
                V.PGG=V.PGG.data, d.V.PGG=d.V.PGG
                
                ) #,



### 2. Initialize model parameters

inits <- function() {
	
#proportional odds model for PGG distribution
rho_int0 <- rnorm((K-1),0,1)
rho_coef <- rnorm(d.V.ETA, mean=0, sd=0.25)
  
#sampling patient's true PGG
eta.hat <- pt.data$bx.pgg[is.na(eta.data)] 

#mixed effects model
xi <- rlnorm(d.Z)
mu_raw <- as.matrix(cbind(rnorm(d.Z),rnorm(d.Z)))
Tau_B_raw <- rwishart((d.Z+1),diag(d.Z)*var_vec)$W
sigma_res <- min(rlnorm(1),1)

#effect of prostate volume
beta <- rnorm(d.X)

#proportional odds model for biopsy results
alpha0 <- rnorm((K-1),0,1)
gamma.PGG <- rnorm((d.V.PGG+(K-1)),mean=0,sd=0.25)
gamma.PGG[(d.V.PGG+1):(d.V.PGG+(K-1))] <- abs(gamma.PGG[(d.V.PGG+1):(d.V.PGG+(K-1))])

list(rho_int0=rho_int0, rho_coef=rho_coef,  
     eta.hat=eta.hat, 
     xi=xi, mu_raw=mu_raw, Tau_B_raw=Tau_B_raw, sigma_res=sigma_res, beta=beta, 
     alpha0=alpha0, gamma.PGG=gamma.PGG
     ) } 


### 3. Define parameters to be tracked
params <- c("rho_int", "rho_coef",
            "eta.hat", 
            "xi", "mu_int", "mu_slope", 
            "sigma_int", "sigma_slope", "sigma_res",  "cov_int_slope", 
            "b.vec", 
            "beta",
            "alpha", "gamma.PGG"
            ) 



