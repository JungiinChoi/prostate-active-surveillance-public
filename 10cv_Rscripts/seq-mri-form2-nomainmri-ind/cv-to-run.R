#Yates Coley
#rycoley@gmail.com
#Updated 2017-7-11
#Annotations updated 12/21/17
#This script will run all scripts for model estimation and preparing patient-level predictions


### WORKFLOW
## 1. Clear workspace
## 2. Define directories, file names
## 3. Load libraries
## 4. Set date of data pull
## 5. Run R scripts for model estimation

### 1. Clear workspace
rm(list=ls())
to.mask<- 1 ## 04/13/22 Zitong

### 2. Define directories, file names
#### These will have to be adjusted by the user
base.location <- '~/Downloads/prostate-active-surveillance-vDaan/' #"/users/rcoley/jhas-epic/psa/psa-new/" #"/Users/ryc/Documents/inhealth/prediction-model/automated/for-TIC/"
location.of.data <- paste0(base.location, "data")
location.of.r.scripts <- paste0(base.location, "10cv_Rscripts/seq-mri-form2-nomainmri")
location.of.generated.files <- paste0(base.location, "generated-files-cv/seq-form2")

name.of.pt.data <- "Demographic_6.15.csv" #"demographics with physician info.2015.csv"
name.of.bx.data <- "Biopsy_6.15.csv" #"Biopsy data_2015.csv"
name.of.psa.data <- "PSA_6.15.csv" #"PSA.2015.csv"
name.of.tx.data <- "Treatment_6.15.csv" #"treatment_2015.csv"

### 3. Load libraries
#### may need to go back and add command for automatic installation
library("lme4")
library("splines")
library("gtools")
library("bayesm")
library("rjags")
library("R2jags")
library("RCurl")


### 4. Set date of data pull
#### This should be adjusted by the user
#date.pull<-as.numeric(as.Date("1/9/18","%m/%d/%y"))
date.pull<-as.numeric(Sys.Date())

### 5. Run R scripts!

#Load data; tidy, check, and shape
#source(paste(location.of.r.scripts,"data-load-check-and-shaping.R",sep="/"))
load("~/Downloads/prostate-active-surveillance-vDaan/generated-files-cv/IOP-data-shaping-work-space.RData")
data.check <- function(condition, message){
  if(condition==FALSE){print(paste(message, "Program terminated.", sep=" "))}
  stopifnot(condition)}

#Load tidied/shaped data; further shaping for JAGS
source(paste(location.of.r.scripts,"cv-data-prep-for-jags-jhpce.R",sep="/"))

#Load arguments for JAGS
source(paste(location.of.r.scripts,"cv-argument-prep-for-jags-jhpce.R",sep="/"))

#Define JAGS model
#source(paste(location.of.r.scripts,"cv-JAGS-prediction-model.R",sep="/"))

#run model and save results
seed <- to.mask
#do.one<-function(seed){
set.seed(seed)
outj<-jags(jags_data, inits=inits, 
           parameters.to.save=params,
           model.file=paste(location.of.r.scripts, 
                            "cv-JAGS-prediction-model-seq-form2.txt", sep="/"),
           n.thin=n.thin, n.chains=n.chains, n.burnin=n.burnin, n.iter=n.iter)
out<-outj$BUGSoutput
  # for(j in 1:length(out$sims.list)){
  # 	write.csv(out$sims.list[[j]],
  # 	          paste(location.of.generated.files, "/jags-prediction-", names(out$sims.list)[j],"-",seed,".csv",sep=""))

write.csv(out$sims.list$eta.track, 
          paste0(location.of.generated.files,"/eta-fitted-",to.mask,".csv"))
  #  }
#}

####this is how I defined the seed on our SGE computing cluster; each task was assigned an id.
# (SEED<-as.numeric(Sys.getenv("SGE_TASK_ID")))
# do.one(seed=SEED)
do.one(2021)

#this will perform model estimation and put posterior estimates in the designated location for generated files

## 6. Run R scripts to process results of model estimation, prep patient predicitons
source(paste(location.of.r.scripts, "patient-predictions-pgs.R", sep="/"))
source(paste(location.of.r.scripts, "get-model-estimation-assessments.R", sep="/"))
source(paste(location.of.r.scripts, "patient-predictions-biopsy-reclassification.R", sep="/"))
source(paste(location.of.r.scripts, "examine-predictions-associations-eta.R", sep="/"))
source(paste(location.of.r.scripts, "examine-predictions-associations.R", sep="/"))

#Import Predictions
source(paste(location.of.r.scripts,"run-prediction-import.R",sep="/"))
