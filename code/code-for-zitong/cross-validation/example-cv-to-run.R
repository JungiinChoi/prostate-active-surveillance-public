# Yates Coley
# rycoley@gmail.com
# 2019 January 08
# This script will run all scripts for model estimation and obtaining predictions for patients with known state


### WORKFLOW
## 1. Clear workspace
## 2. Define directories, file names
## 3. Load libraries
## 4. Set date of data pull
## 5. Set seed, group of patients to mask known state
## 6. Run R scripts to prepare model estimation
## 7. Run model estimation

### 1. Clear workspace
rm(list=ls())

### 2. Define directories, file names
#### These will have to be adjusted by the user 
base.location <- "/users/rcoley/jhas-epic/psa/psa-new/" 
location.of.data <- paste0(base.location, "data")
location.of.r.scripts <- paste0(base.location, "R-scripts")
location.of.generated.files <- paste0(base.location, "generated-files")

name.of.pt.data <- "Demographics_8_9.csv" 
name.of.bx.data <- "Biopsies_8_9.csv" 
name.of.psa.data <- "PSA_8_9.csv" 
name.of.tx.data <- "Treatment_8.9.csv" 



### 3. Load libraries
#### may need to go back and add command for automatic installation
library("lme4")
library("splines")
library("gtools")
library("bayesm")
library("rjags")
library("R2jags")


### 4. Set date of data pull
#### This should be adjusted by the user 
date.pull<-as.numeric(as.Date("8/9/17","%m/%d/%y"))


## 5. Set seed, group of patients to mask known state
#### If using this script to run all R code, can comment out this code for defining the seed and group to be exclude at the top of cv-data-prep-for-jags.R
#### This should be adjusted by the user depending on how seed, masking set 
SEED <- as.numeric(Sys.getenv("SGE_TASK_ID"))
# For GAP3 analysis, make this seed=1,...,10 and define grouping variable with values 1,...,10 
to.mask<- SEED


### 6. Run R scripts to prepare model estimation

#Load data; tidy, check, and shape
source(paste(location.of.r.scripts,"data-load-check-and-shaping.R",sep="/"))

#Load tidied/shaped data; further shaping for JAGS
source(paste(location.of.r.scripts,"cv-data-prep-for-jags.R",sep="/"))

#Load arguments for JAGS
source(paste(location.of.r.scripts,"cv-argument-prep-for-jags.R",sep="/"))

#Define JAGS model
source(paste(location.of.r.scripts,"cv-JAGS-prediction-model.R",sep="/"))

# change length; burn-in; number thinned; number of chains
#### To be adjusted by user 
n.iter <- 40000; n.burnin <- 10000; n.thin <- 20; n.chains <- 1 



## 7. Run model estimation
#run model and save results
set.seed(to.mask)	
outj<-jags(jags_data, inits=inits, 
           parameters.to.save=params, 
           model.file=paste(location.of.r.scripts, "cv-JAGS-prediction-model.txt", sep="/"), 
           n.thin=n.thin, n.chains=n.chains, n.burnin=n.burnin, n.iter=n.iter)

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
