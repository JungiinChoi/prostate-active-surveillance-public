expit <- function(x)
{return(exp(x)/(1+exp(x)))}

library("pROC")
library("ROCR")
library("splines")
location.of.generated.files <- "~/Downloads/prostate-active-surveillance-vDaan/generated-files-seq"
load(paste(location.of.generated.files,"IOP-data-shaping-work-space.RData", sep="/"))

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
my.ci1<-ci.auc(response=RC1, predictor=p_rc1, method="bootstrap") #takes a long time! #0.85, 0.89
my.auc1; my.ci1


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
