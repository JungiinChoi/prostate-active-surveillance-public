#Yates Coley
#rycoley@gmail.com
#Updated 23 August 2018
#This script will run all scripts for model estimation and preparing patient-level predictions


### WORKFLOW
## 1. Clear workspace
## 2. Define directories, file names
## 3. Load libraries
## 4. Set date of data pull
## 5. Run R scripts for model estimation

### 1. Clear workspace
rm(list=ls())

### 2. Define directories, file names
#### These will have to be adjusted by the user 
base.location <- "/users/rcoley/jhas-epic/psa/psa-new/" #"/Users/ryc/Documents/inhealth/prediction-model/automated/for-TIC/"
location.of.data <- paste0(base.location, "data")
location.of.r.scripts <- paste0(base.location, "R-scripts")
location.of.generated.files <- paste0(base.location, "generated-files")

name.of.pt.data <- "Demographics_8_9.csv" #"demographics with physician info.2015.csv"
name.of.bx.data <- "Biopsies_8_9.csv" #"Biopsy data_2015.csv"
name.of.psa.data <- "PSA_8_9.csv" #"PSA.2015.csv"
name.of.tx.data <- "Treatment_8.9.csv" #"treatment_2015.csv"




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

### 5. Run R scripts!

#Load data; tidy, check, and shape
source(paste(location.of.r.scripts,"data-load-check-and-shaping.R",sep="/"))

#Load tidied/shaped data; further shaping for JAGS
source(paste(location.of.r.scripts,"data-prep-for-jags.R",sep="/"))

#Load arguments for JAGS
source(paste(location.of.r.scripts,"argument-prep-for-jags.R",sep="/"))

#Define JAGS model
source(paste(location.of.r.scripts,"JAGS-prediction-model.R",sep="/"))

# change length; burn-in; number thinned; number of chains
n.iter <- 40000; n.burnin <- 10000; n.thin <- 20; n.chains <- 1 

#run model and save results
do.one<-function(seed){
set.seed(seed)	
outj<-jags(jags_data, inits=inits, parameters.to.save=params, 
           model.file=paste(location.of.r.scripts, "JAGS-prediction-model.txt", sep="/"), 
           n.thin=n.thin, n.chains=n.chains, n.burnin=n.burnin, n.iter=n.iter)
out<-outj$BUGSoutput
for(j in 1:length(out$sims.list)){
	write.csv(out$sims.list[[j]], 
	          paste(location.of.generated.files, "/jags-prediction-", names(out$sims.list)[j],"-",seed,".csv",sep=""))} }

(SEED<-as.numeric(Sys.getenv("SGE_TASK_ID"))) 
####this is how I defined the seed on our SGE computing cluster; each task was assigned an id.
do.one(seed=SEED)


#this will perform model estimation and put posterior estimates in the designated location for generated files