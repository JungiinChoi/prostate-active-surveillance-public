## compare sequential model and proportional odds using patients under surgery (eta known)

rm(list=ls())
source('R-scripts/load_libs.R')
base.location <- "/Users/zitongwang/Downloads/prostate-active-surveillance-vDaan/" #"/users/rcoley/jhas-epic/psa/psa-new/" #"/Users/ryc/Documents/inhealth/prediction-model/automated/for-TIC/"
setwd(base.location)
location.of.data <- paste0(base.location, "data")
location.of.r.scripts <- paste0(base.location, "R-scripts")
location.of.generated.files <- paste0(base.location, "generated-files-pgg-model")

name.of.pt.data <- "Demographic_6.15.csv" #"demographics with physician info.2015.csv"
name.of.bx.data <- "Biopsy_6.15.csv" #"Biopsy data_2015.csv"
name.of.psa.data <- "PSA_6.15.csv" #"PSA.2015.csv"
name.of.tx.data <- "Treatment_6.15.csv" #"treatment_2015.csv"

date.pull<-as.numeric(Sys.Date())
source(paste(location.of.r.scripts,"data-load-check-and-shaping.R",sep="/"))
#source(paste(location.of.r.scripts,"data-prep-for-jags-reform.R",sep="/"))
#source(paste(location.of.r.scripts,"argument-prep-for-jags-reform.R",sep="/"))
# source(paste(location.of.r.scripts,"JAGS-prediction-model-etaknown.R",sep="/"))

eta.info <- pt.data %>% dplyr::select(clinical_PTnum,subj, true.pgg)
psa.data.eta <- psa.data %>% left_join(eta.info)
psa.data.eta <- psa.data.eta %>% filter(!is.na(true.pgg))

# use only known eta patients
pt.data.eta1 <- pt.data %>% filter(!is.na(true.pgg))
pt.data.eta2 <- pt.data %>% filter(surg == 1)
pt.data.eta2$id[which(!(pt.data.eta2$id %in% pt.data.eta1$id))] # people went through surgery still has unknown pgg
pt.data.eta <- pt.data.eta1
npat <- dim(pt.data.eta)[1]
cancer_data <- pt.data.eta$true.pgg[!is.na(pt.data.eta$true.pgg)]
npat_cancer_known <- length(cancer_data)

#mean- and varianace- standardized age at diagnosis
pt.data.eta$dx.age.std <- scale(pt.data.eta$dx.age)

modmat_cancer <- as.matrix(cbind(pt.data.eta$dx.age.std, pt.data.eta$lr.vol))
npred_cancer<-dim(modmat_cancer)[2]
#The number of latent classes/ values of true cancer state
nlevel_cancer <- 4
#subset for sequential models
cancer_data1 <- ifelse(cancer_data == 1, 1, 0)
modmat_cancer1 <- modmat_cancer
npat_cancer1 = npat_cancer_known

inx_lev2 <- which(cancer_data > 1)
cancer_data2 <- cancer_data[inx_lev2]
cancer_data2 <- ifelse(cancer_data2 == 2, 1, 0)
modmat_cancer2 <- modmat_cancer[inx_lev2,]
npat_cancer2 = length(cancer_data2)

inx_lev3 <- which(cancer_data > 2)
cancer_data3 <- cancer_data[inx_lev3]
cancer_data3 <- ifelse(cancer_data3 == 3, 1, 0)
modmat_cancer3 <- modmat_cancer[inx_lev3,]
npat_cancer3 = length(cancer_data3)

### marginal probabilities of cancers 
mean(cancer_data1)
mean(cancer_data2)
mean(cancer_data3)

### 1. Set up jags arguments
jags_data<-list(#nlevel_cancer=nlevel_cancer, 
                cancer_data1=cancer_data1, cancer_data2=cancer_data2,
                cancer_data3=cancer_data3,  
                modmat_cancer1=modmat_cancer1, modmat_cancer2=modmat_cancer2,
                modmat_cancer3=modmat_cancer3,  
                npat_cancer1=npat_cancer1,npat_cancer2=npat_cancer2,
                npat_cancer3=npat_cancer3, 
                npred_cancer = npred_cancer,
                alpha = 1, beta = 1
                
) 

### 2. Initialize model parameters

inits <- function() {
  cancer_int1 <- cancer_int2 <-cancer_int3 <-rnorm(1,0,1)
  cancer_slope1 <- rnorm(npred_cancer,mean=0,sd=0.25)
  cancer_slope2 <- rnorm(npred_cancer,mean=0,sd=0.25)
  cancer_slope3 <- rnorm(npred_cancer,mean=0,sd=0.25)
  
  list(
    cancer_int1=cancer_int1, cancer_int2=cancer_int2, cancer_int3=cancer_int3, 
    cancer_slope1=cancer_slope1,cancer_slope2=cancer_slope2,cancer_slope3=cancer_slope3
  ) }


### 3. Define parameters to be tracked
params <- c(
  "cancer_int1", "cancer_int2", "cancer_int3",
  "cancer_slope1", "cancer_slope2", "cancer_slope3")


### 4. Define other jags settings  
n.iter <- 10000; n.burnin <- 2500; n.thin <- 10; n.chains <- 1

### 5. Run
seed=2022
set.seed(seed)
outj<-jags(jags_data, inits=inits, parameters.to.save=params,
           model.file=paste("code/compare_cancerreg_seqJAGS.txt"),
           n.thin=n.thin, n.chains=n.chains, n.burnin=n.burnin, n.iter=n.iter)
out<-outj$BUGSoutput
location.of.generated.files <-  paste0(base.location, "generated-files-pgg-model/sequential_model")
for(j in 1:length(out$sims.list)){
  write.csv(out$sims.list[[j]],
            paste(location.of.generated.files, "/jags-prediction-", names(out$sims.list)[j],"-",seed,".csv",sep=""))
  
}
get_stats <- function(x){a <- quantile(x, c(0.025, 0.975)); b<-mean(x); return(c(a[1], b, a[2]))}


## compare with original code
library(tidyverse)
get_stats <- function(x){a <- quantile(x, c(0.025, 0.975)); b<-mean(x); return(c(a[1], b, a[2]))}
jags_int1 <- 
  read_csv(paste0(location.of.generated.files,"/jags-prediction-cancer_int1-2022.csv"))
jags_int2 <- 
  read_csv(paste0(location.of.generated.files,"/jags-prediction-cancer_int2-2022.csv"))
jags_int3 <- 
  read_csv(paste0(location.of.generated.files,"/jags-prediction-cancer_int3-2022.csv"))
apply(jags_int1, 2, get_stats)
apply(jags_int2, 2, get_stats)
apply(jags_int3, 2, get_stats)

jags_slope1 <- 
  read_csv(paste0(location.of.generated.files,"/jags-prediction-cancer_slope1-2022.csv"))
jags_slope2 <- 
  read_csv(paste0(location.of.generated.files,"/jags-prediction-cancer_slope2-2022.csv"))
jags_slope3 <- 
  read_csv(paste0(location.of.generated.files,"/jags-prediction-cancer_slope3-2022.csv"))
apply(jags_slope1, 2, get_stats)
apply(jags_slope2, 2, get_stats)
apply(jags_slope3, 2, get_stats)

## check calibration
jags_int1 <- as.matrix(jags_int1[,-1])
jags_int2 <- as.matrix(jags_int2[,-1])
jags_int3 <- as.matrix(jags_int3[,-1])
jags_slope1 <- as.matrix(jags_slope1[,-1])
jags_slope2 <- as.matrix(jags_slope2[,-1])
jags_slope3 <- as.matrix(jags_slope3[,-1])
expit <- function(x){exp(x)/(1+exp(x))}

nsim <- ((n.iter - n.burnin)/n.thin)
p_rc1 <- p_rc2 <- p_rc3 <-  p_rc4 <- matrix(0, ncol = nsim, nrow = npat_cancer_known)
for(i in 1:nsim){
  int1 <- matrix(jags_int1[i,])
  int2 <-  matrix(jags_int2[i,])
  int3 <-  matrix(jags_int3[i,])
  slope1 <- matrix(jags_slope1[i,])
  slope2 <- matrix(jags_slope2[i,])
  slope3 <- matrix(jags_slope3[i,])
  
  #prediction using all data (modmat_cancer)
  linpred1 <- modmat_cancer %*% slope1
  p_y1_given_x0 <- expit(matrix(int1 + c(linpred1)))
  linpred2 <- modmat_cancer %*% slope2
  p_y2_given_x0ge1 <- expit(matrix(int2 + c(linpred2)))
  linpred3 <- modmat_cancer %*% slope3
  p_y3_given_x0ge2 <- expit(matrix(int3 + c(linpred3)))
  
  p_rc1[,i] <- p_y1_given_x0
  p_rc2[,i] <- (1-p_y1_given_x0) * p_y2_given_x0ge1
  p_rc3[,i] <- (1-p_y1_given_x0)*(1-p_y2_given_x0ge1) * p_y3_given_x0ge2
  p_rc4[,i] <- 1- p_rc1[,i] -  p_rc2[,i] -  p_rc3[,i]
}

obs1 <- sum(cancer_data == 1);obs2 <- sum(cancer_data == 2)
obs3 <- sum(cancer_data == 3);obs4 <- sum(cancer_data == 4)

calib <-  matrix(0, nrow=nsim)
exp1 <- exp2 <- exp3 <-exp4 <- calib <- matrix(0, nrow = nsim)
# for(i in 1:nsim){
#   exp1[i] <- sum(p_rc1[,i]);exp2[i] <- sum(p_rc2[,i])
#   exp3[i] <- sum(p_rc3[,i]);exp4[i] <- sum(p_rc4[,i])
# 
#   calib[i] <-
#     (obs1-exp1[i])^2/exp1[i] + (obs2-exp2[i])^2/exp2[i] +
#     (obs3-exp3[i])^2/exp3[i] + (obs4-exp4[i])^2/exp4[i]
# }
# quantile(calib, c(0.025, 0.5, 0.975))
# summary(calib)


exp1 <- sum(apply(p_rc1, 1, mean))
exp2 <- sum(apply(p_rc2, 1, mean))
exp3 <- sum(apply(p_rc3, 1, mean))
exp4 <- sum(apply(p_rc4, 1, mean))

(obs1-exp1)^2/exp1 + (obs2-exp2)^2/exp2 + 
  (obs3-exp3)^2/exp3 + (obs4-exp4)^2/exp4

