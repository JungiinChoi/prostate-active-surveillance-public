#Yates Coley
#rycoley@gmail.com
#2019 January 08
#This script runs the JAGS model for 10-fold CV

#### WORKFLOW
### 0. Load R packages
### 1. Define JAGS model
### 2. Define estimation settings ### *Daan* These will need to be changed.
### 3. Run JAGS estimation and save results ### *Daan* Randomization starting seed will need to be changed


### 0. Load R packages
library("gtools") #use logit() from this package
library("rjags")
library("R2jags")



### 1. Define JAGS model
source(paste(location.of.r.scripts,"cv-JAGS-prediction-model.R",sep="/"))



### 2. Define estimation settings ### *Daan* These will need to be changed.

# The first time you run this code, deine a shorter number of iterations to check that the code compiles
n.iter <- 100; n.burnin <- 20; n.thin <- 5; n.chains <- 1 

# After code has successfull compiled, you can comment out shorter run above and set a longer run
#n.iter <- 10000; n.burnin <- 1000; n.thin <- 20; n.chains <- 1 

#I genuniely have no idea how long this model will need to run. 
#There are auto-update options (autojags()) we can use to run the posterior sample if it hasn't converged.

#Note that I only run 1 chain here. We want to run 5 chains total. In my work, I run them in parallel. 



### 3. Run JAGS estimation and save results ### *Daan* Seed will need to be changed

#SEED was previously defined here, it is now defined at the top of cv-data-prep-for-jags.R
#SEED=to.mask
set.seed(SEED)	


# run JAGS model
outj<-jags(jags_data, inits=inits, parameters.to.save=params, 
           model.file=paste(location.of.r.scripts, "cv-JAGS-prediction-model.txt", sep="/"), 
           n.thin=n.thin, n.chains=n.chains, n.burnin=n.burnin, n.iter=n.iter)

#update until convergence
#  outj <- autojags(outj) #we can use this if the posterior is not converging and we need to run the chain longer

#save results
out<-outj$BUGSoutput
write.csv(out$sims.list$eta.track, 
          paste0(location.of.generated.files,"/eta-fitted-",to.mask,".csv"))


#This script saves additional parameters (beyond patient-specific eta predictions) 
#in order to monitor convergence across many chains 

for(j in 1:length(out$sims.list)){
  write.csv(out$sims.list[[j]], 
            paste(location.of.generated.files, 
                  "/cv-jags-prediction-", 
                  names(out$sims.list)[j],"-", 
                  to.mask,
                  ".csv", sep=""))} 
