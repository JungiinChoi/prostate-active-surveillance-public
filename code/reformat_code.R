rm(list=ls())
setwd("~/Downloads/prostate-active-surveillance-vDaan")
source('R-scripts/load_libs.R')
location.of.data <- "data"
location.of.r.scripts <- "R-scripts"
location.of.generated.files <- "generated-files-small"

name.of.pt.data <- "Demographic_6.15.csv" #"demographics with physician info.2015.csv"
name.of.bx.data <- "Biopsy_6.15.csv" #"Biopsy data_2015.csv"
name.of.psa.data <- "PSA_6.15.csv" #"PSA.2015.csv"
name.of.tx.data <- "Treatment_6.15.csv" #"treatment_2015.csv"

pt.data.original <- read.csv(paste0(location.of.data, "/",name.of.pt.data))
bx.data.original <- read.csv(paste0(location.of.data, "/",name.of.bx.data))
psa.data.original <- read.csv(paste0(location.of.data, "/",name.of.psa.data))
tx.data.original <- read.csv(paste0(location.of.data, "/",name.of.tx.data))

## subset 100 patient (use psa.data patients) # which(!(psa.data.original$clinical_PTnum %in% pt.data.original$clinical_PTnum))
nsub <- 150
ptnum <- unique(psa.data.original$clinical_PTnum)
set.seed(2021)
sub_ptnum <- sample(ptnum, nsub, replace = F)
pt.data <- pt.data.original[pt.data.original$clinical_PTnum %in% sub_ptnum,]

## match other datasets
bx.data <-  bx.data.original[bx.data.original$clinical_PTnum %in% sub_ptnum,]
psa.data <-  psa.data.original[psa.data.original$clinical_PTnum %in% sub_ptnum,]
tx.data <-  tx.data.original[tx.data.original$clinical_PTnum %in% sub_ptnum,]


write.csv(psa.data, "generated-files-small/PSA_subset.csv", row.names = F)
write.csv(pt.data, "generated-files-small/Demographic_subset.csv", row.names = F)
write.csv(bx.data, "generated-files-small/Biopsy_subset.csv", row.names = F)
write.csv(tx.data, "generated-files-small/Treatment_subset.csv", row.names = F)


date.pull<-as.numeric(Sys.Date())
name.of.pt.data <- "Demographic_subset.csv" #"demographics with physician info.2015.csv"
name.of.bx.data <- "Biopsy_subset.csv" #"Biopsy data_2015.csv"
name.of.psa.data <- "PSA_subset.csv" #"PSA.2015.csv"
name.of.tx.data <- "Treatment_subset.csv" #"treatment_2015.csv"

location.of.data <- "generated-files-small"
location.of.generated.files <- "generated-files-small/with_Yates_code/"

source(paste(location.of.r.scripts,"data-load-check-and-shaping.R",sep="/"))
source(paste(location.of.r.scripts,"data-prep-for-jags-reform.R",sep="/"))
source(paste(location.of.r.scripts,"argument-prep-for-jags-reform.R",sep="/"))
source(paste(location.of.r.scripts,"JAGS-prediction-model-reform.R",sep="/"))


seed <- 2021
#do.one<-function(seed){
  set.seed(seed)
  outj<-jags(jags_data, inits=inits, parameters.to.save=params,
             model.file=paste(location.of.r.scripts, "JAGS-prediction-model-reform.txt", sep="/"),
             n.thin=n.thin, n.chains=n.chains, n.burnin=n.burnin, n.iter=n.iter)
  out<-outj$BUGSoutput
  for(j in 1:length(out$sims.list)){
    write.csv(out$sims.list[[j]],
              paste(location.of.generated.files, "/with_reformat_code/jags-prediction-", names(out$sims.list)[j],"-",seed,".csv",sep=""))
  }
#}


#do.one(2021)

library(readr)
xi_yates <- read_csv("generated-files-small/with_Yates_code/jags-prediction-xi-2021.csv")
xi_reform <- read_csv("generated-files-small/with_reformat_code/jags-prediction-xi-2021.csv")
head(xi_yates)
head(xi_reform)
