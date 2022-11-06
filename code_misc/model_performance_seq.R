expit <- function(x)
{return(exp(x)/(1+exp(x)))}

library("pROC")
library("ROCR")
library("splines")
library("lme4")
library("gtools")
library("bayesm")
library("rjags")
library("R2jags")
library("RCurl")
library("readr")
library("dplyr")
location.of.r.scripts <- paste0("~/Downloads/prostate-active-surveillance-vDaan/R-script")
location.of.generated.files <- "~/Downloads/prostate-active-surveillance-vDaan/generated-files-seq-mri-form2"
load(paste(location.of.generated.files,"IOP-data-shaping-work-space.RData", sep="/"))
options(warn=1)
data.check <- function(condition, message){
  if(condition==FALSE){print(paste(message, "Program terminated.", sep=" "))}
  stopifnot(condition)}
source(paste(location.of.r.scripts,"data-prep-for-jags-seq-mri-form2.R",sep="/"))
#number of patients
(n <- dim(pt.data)[1])

#get observed PGGs
eta.data<-pt.data$true.pgg
obs.eta<-eta.data[!is.na(eta.data)]
#put in matrix for PGG>1, PGG>2, PGG>3
obs.eta.mat<-cbind(as.numeric(obs.eta>1), as.numeric(obs.eta>2), as.numeric(obs.eta>3))

#etahat<-read.csv(paste0(location.of.generated.files, "/jags-prediction-eta.hat-2021.csv"))
etahat<-read.csv(paste0(location.of.generated.files, "/jags-prediction-cancer_state-2022.csv"))
etahat<-as.matrix(etahat[,2:dim(etahat)[2]])

K<-4 #number of classes
(B<-dim(etahat)[1]) #length of saved posterior chain
(N<-length(obs.eta))
npat_cancer_known <- length(obs.eta)
#mean eta predictions for patients with unknown cancer state
pred_eta <- matrix(nrow=(n-npat_cancer_known), ncol=K)
for(i in 1:(n-npat_cancer_known)){for(k in 1:K){
  pred_eta[i,k] <- sum(etahat[,i]==k)/B }}

## get parameter estimates in biopsy outcome proportional odds model
#alpha<-read.csv(paste0(location.of.generated.files, "/jags-prediction-alpha-2021.csv"))
alpha1<-read.csv(paste0(location.of.generated.files, "/jags-prediction-pgg_int1-2022.csv"))
alpha1<-as.matrix(alpha1[,2:dim(alpha1)[2]])
# posterior estimate
alpha1_med<-apply(alpha1,2,median)

alpha2<-read.csv(paste0(location.of.generated.files, "/jags-prediction-pgg_int2-2022.csv"))
alpha2<-as.matrix(alpha2[,2:dim(alpha2)[2]])
# posterior estimate
alpha2_med<-apply(alpha2,2,median)

alpha3<-read.csv(paste0(location.of.generated.files, "/jags-prediction-pgg_int3-2022.csv"))
alpha3<-as.matrix(alpha3[,2:dim(alpha3)[2]])
# posterior estimate
alpha3_med<-apply(alpha3,2,median)

alpha_med <- c(alpha1_med, alpha2_med, alpha3_med)

gamma1<-read.csv(paste0(location.of.generated.files, "/jags-prediction-pgg_slope1-2022.csv"))
gamma1<-as.matrix(gamma1[,2:dim(gamma1)[2]])
# posterior estimate
gamma1_med<-apply(gamma1,2,median)

gamma2<-read.csv(paste0(location.of.generated.files, "/jags-prediction-pgg_slope2-2022.csv"))
gamma2<-as.matrix(gamma2[,2:dim(gamma2)[2]])
# posterior estimate
gamma2_med<-apply(gamma2,2,median)

gamma3<-read.csv(paste0(location.of.generated.files, "/jags-prediction-pgg_slope3-2022.csv"))
gamma3<-as.matrix(gamma3[,2:dim(gamma3)[2]])
# posterior estimate
gamma3_med<-apply(gamma3,2,median)

gamma_med <- rbind(gamma1_med, gamma2_med, gamma3_med)

## update covariate matrix for biopsy model with eta predictions
V.PGG.eta<-matrix(nrow=npat_pgg, ncol=(K-1))
for(j in 1:npat_pgg){
  if(pgg_patient_index_map[j] <= npat_cancer_known){
    for(k in 1:(K-1)){
      V.PGG.eta[j,k] <- as.numeric(eta.data[pgg_patient_index_map[j]]>k)}}
  else{
    for(k in 1:(K-1)){
      V.PGG.eta[j,k] <- 1-sum(pred_eta[(pgg_patient_index_map[j]-npat_cancer_known),1:k]) } } }

## get estimated probability of upgrading
## for PGG >1 ------------------
#P(pgg = 1|x)
lin_pred1 <- alpha_med[1] + as.matrix(cbind(modmat_pgg, V.PGG.eta)) %*% as.vector(gamma_med[1,])
#on probability scale
p_rc1 <- 1-expit(lin_pred1)

#observed outcome
RC1 <- as.numeric(pgg_data>1)

#AUC for RC predictions
my.auc1<-performance(prediction(p_rc1, RC1), "auc")@y.values[[1]]  #0.887
#my.ci1<-ci.auc(response=RC1, predictor=p_rc1, method="bootstrap") #takes a long time! #0.85, 0.89
my.auc1 


## for PGG >2 ------------------
#linear predictor
lin_pred2 <- alpha_med[2] + as.matrix(cbind(modmat_pgg, V.PGG.eta)) %*% as.vector(gamma_med[2,])
#on probability scale
p_rc2 <- 1- 1/(1+exp(lin_pred1)) * expit(lin_pred2) - expit(lin_pred1)

#observed outcome
RC2 <- as.numeric(pgg_data>2)

#AUC for RC predictions
my.auc2<-performance(prediction(p_rc2, RC2), "auc")@y.values[[1]]  #0.887
my.ci2<-ci.auc(response=RC2, predictor=p_rc2, method="bootstrap") #takes a long time! #0.85, 0.89
my.auc2; my.ci2

## for PGG >3 ------------------
#linear predictor
lin_pred3 <- alpha_med[3] + as.matrix(cbind(modmat_pgg, V.PGG.eta)) %*% as.vector(gamma_med[3,])
#on probability scale
p_rc3 <- 1-1/(1+exp(lin_pred1))*1/(1+exp(lin_pred2))*expit(lin_pred3) - 1/(1+exp(lin_pred1)) * expit(lin_pred2) - expit(lin_pred1)

#observed outcome
RC3 <- as.numeric(pgg_data>3)

#AUC for RC predictions
my.auc3<-performance(prediction(p_rc3, RC3), "auc")@y.values[[1]]  #0.887
my.ci3<-ci.auc(response=RC3, predictor=p_rc3, method="bootstrap") #takes a long time! #0.85, 0.89
my.auc3;my.ci3


### Calibration for PGG
exp1_pgg <- sum(1-p_rc1)
exp2_pgg <- sum(p_rc1 -p_rc2)
exp3_pgg <- sum(p_rc2 -p_rc3)
exp4_pgg <- sum(p_rc3)
get_chisq <- function(obs, exp){
  (obs[1]-exp[1])^2/exp[1] + (obs[2]-exp[2])^2/exp[2] + 
    (obs[3]-exp[3])^2/exp[3] + (obs[4]-exp[4])^2/exp[4]
}
obs1_pgg <- sum(pgg_data == 1)
obs2_pgg <- sum(pgg_data == 2)
obs3_pgg <- sum(pgg_data == 3)
obs4_pgg <- sum(pgg_data == 4)
obs_pgg <- c(obs1_pgg,obs2_pgg,obs3_pgg, obs4_pgg)
exp_pgg <- c(exp1_pgg,exp2_pgg,exp3_pgg, exp4_pgg)
get_chisq(obs_pgg, exp_pgg)
#### Check eta ------------
id_eta_obs <- pt.data$clinical_PTnum[!is.na(pt.data$true.pgg)]
eta_obs <- pt.data$true.pgg[!is.na(pt.data$true.pgg)]
inx_eta_obs <- which(!is.na(pt.data$true.pgg))
V.ETA.data.obs <- V.ETA.data[inx_eta_obs,]
#### sequential replaced model
cancer_int1 <- read_csv("generated-files-seq-addmri/jags-prediction-cancer_int1-2022.csv")
cancer_int2 <- read_csv("generated-files-seq-addmri/jags-prediction-cancer_int2-2022.csv")
cancer_int3 <- read_csv("generated-files-seq-addmri/jags-prediction-cancer_int3-2022.csv")

cancer_slope1 <- read_csv("generated-files-seq-addmri/jags-prediction-cancer_slope1-2022.csv")
cancer_slope2 <- read_csv("generated-files-seq-addmri/jags-prediction-cancer_slope2-2022.csv")
cancer_slope3 <- read_csv("generated-files-seq-addmri/jags-prediction-cancer_slope3-2022.csv")

cancer_int1 <- cancer_int1[,-1];cancer_int2 <- cancer_int2[,-1];cancer_int3 <- cancer_int3[,-1]
cancer_slope1 <- cancer_slope1[,-1];cancer_slope2 <- cancer_slope2[,-1];cancer_slope3 <- cancer_slope3[,-1]

int1_mean <-apply(cancer_int1, 2, mean) 
int2_mean <-apply(cancer_int2, 2, mean) 
int3_mean <-apply(cancer_int3, 2, mean) 

slope1_mean <-apply(cancer_slope1, 2, mean) 
slope2_mean <-apply(cancer_slope2, 2, mean) 
slope3_mean <-apply(cancer_slope3, 2, mean) 

linpred1 <- int1_mean +  V.ETA.data.obs %*% slope1_mean
p1_eta_seq <- exp(linpred1)/(1+exp(linpred1))

linpred2 <- int2_mean +  V.ETA.data.obs %*% slope2_mean
cond_p2_eta <- exp(linpred2)/(1+exp(linpred2))
p2_eta_seq <- 1/(1+exp(linpred1)) * cond_p2_eta

linpred3 <- int3_mean +  V.ETA.data.obs %*% slope3_mean
cond_p3_eta <- exp(linpred3)/(1+exp(linpred3))
p3_eta_seq <- 1/(1+exp(linpred1)) * 1/(1+exp(linpred2)) * cond_p3_eta

p4_eta_seq <- 1-p1_eta_seq - p2_eta_seq - p3_eta_seq
exp1 <- sum(p1_eta_seq)
exp2 <- sum(p2_eta_seq)
exp3 <- sum(p3_eta_seq)
exp4 <- sum(p4_eta_seq)
exp_eta_seq <- c(exp1, exp2, exp3, exp4)

obs1 <- sum(eta_obs==1)
obs2 <- sum(eta_obs==2)
obs3 <- sum(eta_obs==3)
obs4 <- sum(eta_obs==4)
obs_eta_original <- c(obs1, obs2, obs3, obs4)
get_chisq(obs_eta_original, exp_eta_seq)

obs1_eta <- as.numeric(eta_obs <= 2)
auc(obs1_eta, p1_eta_seq + p2_eta_seq)


### How cancer affect pirads
pirads_int <- read_csv("generated-files-seq-addmri/jags-prediction-pirads_int-2022.csv")
pirads_slope <- read_csv("generated-files-seq-addmri/jags-prediction-pirads_slope-2022.csv")
apply(pirads_int,2, mean)
apply(pirads_slope,2, mean)



