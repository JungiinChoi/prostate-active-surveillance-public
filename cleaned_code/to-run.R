rm(list=ls())

mri_role <- "both" #moderator, both, outcome, 0

# base.location <- "/users/zwang3/PAS/"
# location.of.data <- paste0(base.location, "data")
# location.of.r.scripts <- paste0(base.location, "code")
# location.of.generated.files <- paste0(base.location, "generated-files")

base.location <- "/Users/zitongwang/Downloads/prostate-active-surveillance-vDaan/"
location.of.data <- paste0(base.location, "data")
location.of.r.scripts <- paste0(base.location, "cleaned_code")
location.of.generated.files <- paste0(base.location, "cleaned_code/generated-files")


name.of.pt.data <- "Demographic_6.15.csv" #"demographics with physician info.2015.csv"
name.of.bx.data <- "Biopsy_6.15.csv" #"Biopsy data_2015.csv"
name.of.psa.data <- "PSA_6.15.csv" #"PSA.2015.csv"
name.of.tx.data <- "Treatment_6.15.csv" #"treatment_2015.csv"
name.of.mri.data <- "MRI_6.15.csv" 

library("lme4")
library("splines")
library("gtools")
library("bayesm")
library("rjags")
library("R2jags")
library("RCurl")
library("tidyverse")
library("dplyr")
library("readr")

### 5. Run R scripts!
data.check <- function(condition, message){
  if(condition==FALSE){print(paste(message, "Program terminated.", sep=" "))}
  stopifnot(condition)}

#Load data; tidy, check, and shape
load(paste(location.of.r.scripts,"IOP-data-shaping-work-space.RData",sep="/"))

#Load tidied/shaped data; further shaping for JAGS
source(paste(location.of.r.scripts,"data-prep-for-jags.R",sep="/"))

#Load arguments for JAGS
source(paste(location.of.r.scripts,"argument-prep-for-jags.R",sep="/"))

#Define JAGS model
source(paste(location.of.r.scripts,"JAGS-prediction-model.R",sep="/"))

#run model and save results
seed <- 2022
#do.one<-function(seed){
set.seed(seed)
outj<-jags(jags_data, inits=inits, parameters.to.save=params,
           model.file=paste(location.of.r.scripts, 
                            paste0("JAGS-prediction-model-",mri_role,".txt"), sep="/"),
           n.thin=n.thin, n.chains=n.chains, n.burnin=n.burnin, n.iter=n.iter)
out<-outj$BUGSoutput
for(j in 1:length(out$sims.list)){
  write.csv(out$sims.list[[j]],
            paste(location.of.generated.files, 
                  "/jags-prediction-", 
                  names(out$sims.list)[j],"-",seed,"-", 
                  mri_role, ".csv",sep=""))
}
#}



