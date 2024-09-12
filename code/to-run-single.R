#Yates Coley
#rycoley@gmail.com
#This script will run all scripts for model estimation and preparing patient-level predictions


### WORKFLOW
## 1. Clear workspace
## 2. Define directories, file names
## 3. Load libraries
## 4. Set date of data pull
## 5. Run R scripts for model estimation

### 1. Clear workspace
rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)
mri_role <- args[1]
workdir <- args[2]
J <- args[3]

### 2. Define directories, file names

base.location <- workdir
location.of.data <- paste0(base.location, "/data")
location.of.r.scripts <- paste0(base.location, "/code")
location.of.generated.files <- paste0(base.location, "/generated-files")
location.of.generated.folder = paste(location.of.generated.files, "/", mri_role, sep="")
ifelse(!dir.exists(location.of.generated.folder), dir.create(location.of.generated.folder), FALSE)

### 3. Load libraries
#### may need to go back and add command for automatic installation 
# Package names
packages <- c("lme4",  "dplyr", "tidyr", "readr",
              "splines", "bayesm", "rjags", "R2jags")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

### 5. Run R scripts!

#Load data; tidy, check, and shape
data.check <- function(condition, message){
  if(condition==FALSE){print(paste(message, "Program terminated.", sep=" "))}
  stopifnot(condition)}

#Load data; tidy, check, and shape
options(warn = 0)
#Load tidied/shaped data; further shaping for JAGS
source(paste(location.of.r.scripts,"data-prep-for-jags.R",sep="/"))

#Load arguments for JAGS
source(paste(location.of.r.scripts,"argument-prep-for-jags.R",sep="/"))

#run model and save results
seed <- 2024

#Define JAGS model
source(paste(location.of.r.scripts,"JAGS-prediction-model.R",sep="/"))

set.seed(seed)

outj<-R2jags::jags(jags_data, inits=inits, 
           parameters.to.save=params,
           model.file=paste(location.of.r.scripts, 
                            jags_file_name, sep="/"),
           n.thin=n.thin, n.chains=n.chains, n.burnin=n.burnin, n.iter=n.iter)
out<-outj$BUGSoutput


for(i in 1:length(out$sims.list)){
  if (names(out$sims.list)[i] %in% "p_eta") {
    tmp_mean <- apply(out$sims.list[[i]], 2:4, mean)
    for (k in 1:J){
      write.csv(tmp_mean[k,,],
                paste(location.of.generated.folder, "/jags-prediction-", names(out$sims.list)[i],"_", k, "-", mri_role,".csv",sep=""))
    }
  }
}

