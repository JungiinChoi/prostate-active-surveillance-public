# Yates Coley
# rycoley@gmail.com
# 2019 January 07
# This script prepares model summaries and assessments for the primary model estimation.



#### WORKFLOW
### 1. Clear workspace
### 2. Define directories, file names #### User will need to edit
### 3. Load libraries
### 4. Set date of data pull #### User will need to edit
### 5. Load data
### 6. Run R scripts!


### 1. Clear workspace
rm(list=ls())


### 2. Define directories, file names
#### User will need to edit
base.location <- "/Users/zitongwang/Downloads/prostate-active-surveillance-vDaan/" #"/users/rcoley/jhas-epic/psa/psa-new/" #"/Users/ryc/Documents/inhealth/prediction-model/automated/for-TIC/"
location.of.data <- paste0(base.location, "data")
location.of.r.scripts <- paste0(base.location, "R-scripts")
location.of.generated.files <- paste0(base.location, "generated-files")

#consider separate location for assessment summary R scripts
location.of.assessment.r.scripts<- paste0(base.location, "code/code-for-zitong/model-diagnostics")
location.of.assessment.summaries <- location.of.generated.files


#### User will need to edit
name.of.pt.data <- "Demographics_6.15.csv"  #"Demographics_8_9.csv" 
name.of.bx.data <- "Biopsies_6.15.csv"
name.of.psa.data <- "PSA_6.15.csv"
name.of.tx.data <- "Treatment_6.15.csv"


  
### 3. Load libraries
library("mixtools")
library("lme4")
library("splines")
library("gtools")
library("bayesm")
library("coda")
library("ROCR")
library("pROC")



### 4. Set date of data pull
#### User will need to edit
date.pull<-as.numeric(as.Date("03/21/22","%m/%d/%y"))



### 5. Load data
#Load data; tidy, check, and shape
#source(paste(location.of.r.scripts,"data-load-check-and-shaping.R",sep="/"))

#it is faster to load RData workspace than data-load-check-and-shaping.R
#define data.check function
data.check <- function(condition, message){
  if(condition==FALSE){print(paste(message, "Program terminated.", sep=" "))}
  stopifnot(condition)}

load(paste0(location.of.generated.files,"/IOP-data-shaping-work-space.RData"))


#Load tidied/shaped data; further shaping for JAGS
source(paste(location.of.r.scripts,"data-prep-for-jags.R",sep="/"))

#Load arguments for JAGS
source(paste(location.of.r.scripts,"argument-prep-for-jags.R",sep="/"))

#Define JAGS model
source(paste(location.of.r.scripts,"JAGS-prediction-model.R",sep="/"))




### 6. Run R scripts for assessment

# calibration plot for observables
source(paste0(location.of.assessment.r.scripts,"/calibration-plot-observables.R"))

# convergence diagnostics
source(paste0(location.of.assessment.r.scripts,"/convergence-diagnostics.R"))

# PSA random effects estimates
source(paste0(location.of.assessment.r.scripts,"/psa-lme-plot.R"))

