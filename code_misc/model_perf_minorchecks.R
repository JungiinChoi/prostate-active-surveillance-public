expit <- function(x)
{return(exp(x)/(1+exp(x)))}

library("pROC")
library("ROCR")
library("splines")
location.of.generated.files <- "~/Downloads/prostate-active-surveillance-vDaan/generated-files-seq"
load(paste(location.of.generated.files,"IOP-data-shaping-work-space.RData", sep="/"))
source(paste(location.of.r.scripts,"data-prep-for-jags-reform.R",sep="/"))

#number of patients
(n <- dim(pt.data)[1])

#get observed PGGs
eta.data<-pt.data$true.pgg
obs.eta<-eta.data[!is.na(eta.data)]
#put in matrix for PGG>1, PGG>2, PGG>3
obs.eta.mat<-cbind(as.numeric(obs.eta>1), as.numeric(obs.eta>2), as.numeric(obs.eta>3))

#etahat<-read.csv(paste0(location.of.generated.files, "/jags-prediction-eta.hat-2021.csv"))
etahat_seq<-read.csv(paste0(location.of.generated.files, "/jags-prediction-cancer_state-2022.csv"))
etahat_seq<-as.matrix(etahat_seq[,2:dim(etahat_seq)[2]])
#mean eta predictions for patients with unknown cancer state
pred_eta_seq <- matrix(nrow=(n-npat_cancer_known), ncol=K)
for(i in 1:(n-npat_cancer_known)){for(k in 1:K){
  pred_eta_seq[i,k] <- sum(etahat_seq[,i]==k)/B }}

etahat_canc<-read.csv(paste0("~/Downloads/prostate-active-surveillance-vDaan/generated-files-cancerpom-replaced", "/jags-prediction-cancer_state-2022.csv"))
etahat_canc<-as.matrix(etahat_canc[,2:dim(etahat_canc)[2]])
#mean eta predictions for patients with unknown cancer state
pred_eta_canc <- matrix(nrow=(n-npat_cancer_known), ncol=K)
for(i in 1:(n-npat_cancer_known)){for(k in 1:K){
  pred_eta_canc[i,k] <- sum(etahat_canc[,i]==k)/B }}

etahat<-read.csv(paste0("~/Downloads/prostate-active-surveillance-vDaan/generated-files", "/jags-prediction-eta.hat-2021.csv"))
etahat<-as.matrix(etahat[,2:dim(etahat)[2]])
#mean eta predictions for patients with unknown cancer state
pred_eta <- matrix(nrow=(n-npat_cancer_known), ncol=K)
for(i in 1:(n-npat_cancer_known)){for(k in 1:K){
  pred_eta[i,k] <- sum(etahat[,i]==k)/B }}


apply(pred_eta, 2, sum)
apply(pred_eta_canc, 2, sum)
apply(pred_eta_seq, 2, sum)


## check PGG values for seq model-----------------
location.of.generated.files <- "~/Downloads/prostate-active-surveillance-vDaan/generated-files-seq"
load(paste(location.of.generated.files,"IOP-data-shaping-work-space.RData", sep="/"))
(n <- dim(pt.data)[1])
#get observed PGGs
eta.data<-pt.data$true.pgg
obs.eta<-eta.data[!is.na(eta.data)]
obs.eta.mat<-cbind(as.numeric(obs.eta==1), as.numeric(obs.eta==2), as.numeric(obs.eta==3))
etahat<-read.csv(paste0(location.of.generated.files, "/jags-prediction-cancer_state-2022.csv"))
etahat<-as.matrix(etahat[,2:dim(etahat)[2]])
K<-4 
(B<-dim(etahat)[1]) 
(N<-length(obs.eta))
pred_eta <- matrix(nrow=(n-npat_cancer_known), ncol=K)
for(i in 1:(n-npat_cancer_known)){for(k in 1:K){
  pred_eta[i,k] <- sum(etahat[,i]==k)/B }}
alpha1<-read.csv(paste0(location.of.generated.files, "/jags-prediction-pgg_int1-2022.csv"))
alpha1<-as.matrix(alpha1[,2:dim(alpha1)[2]])
alpha1_med<-apply(alpha1,2,median)
alpha2<-read.csv(paste0(location.of.generated.files, "/jags-prediction-pgg_int2-2022.csv"))
alpha2<-as.matrix(alpha2[,2:dim(alpha2)[2]])
alpha2_med<-apply(alpha2,2,median)
alpha3<-read.csv(paste0(location.of.generated.files, "/jags-prediction-pgg_int3-2022.csv"))
alpha3<-as.matrix(alpha3[,2:dim(alpha3)[2]])
alpha3_med<-apply(alpha3,2,median)
alpha_med <- c(alpha1_med, alpha2_med, alpha3_med)
gamma1<-read.csv(paste0(location.of.generated.files, "/jags-prediction-pgg_slope1-2022.csv"))
gamma1<-as.matrix(gamma1[,2:dim(gamma1)[2]])
gamma1_med<-apply(gamma1,2,median)
gamma2<-read.csv(paste0(location.of.generated.files, "/jags-prediction-pgg_slope2-2022.csv"))
gamma2<-as.matrix(gamma2[,2:dim(gamma2)[2]])
gamma2_med<-apply(gamma2,2,median)
gamma3<-read.csv(paste0(location.of.generated.files, "/jags-prediction-pgg_slope3-2022.csv"))
gamma3<-as.matrix(gamma3[,2:dim(gamma3)[2]])
gamma3_med<-apply(gamma3,2,median)
gamma_med <- rbind(gamma1_med, gamma2_med, gamma3_med)
V.PGG.eta<-matrix(nrow=npat_pgg, ncol=(K-1))
for(j in 1:npat_pgg){
  if(pgg_patient_index_map[j] <= npat_cancer_known){
    for(k in 1:(K-1)){
      V.PGG.eta[j,k] <- as.numeric(eta.data[pgg_patient_index_map[j]]>k)}}
  else{
    for(k in 1:(K-1)){
      V.PGG.eta[j,k] <- 1-sum(pred_eta[(pgg_patient_index_map[j]-npat_cancer_known),1:k]) } } }
lin_pred1 <- alpha_med[1] + as.matrix(cbind(modmat_pgg, V.PGG.eta)) %*% as.vector(gamma_med[1,])
lin_pred2 <- alpha_med[2] + as.matrix(cbind(modmat_pgg, V.PGG.eta)) %*% as.vector(gamma_med[2,])
lin_pred3 <- alpha_med[3] + as.matrix(cbind(modmat_pgg, V.PGG.eta)) %*% as.vector(gamma_med[3,])
exp1_seq <- expit(lin_pred1)
exp2_seq <- (1-expit(lin_pred1))*expit(lin_pred2)
exp3_seq <- (1-expit(lin_pred1))*(1-expit(lin_pred2))*expit(lin_pred3)
exp4_seq <- 1-exp1_seq-exp2_seq-exp3_seq
exp_seq <- cbind(exp1_seq,exp2_seq,exp3_seq,exp4_seq)
apply(exp_seq, 2, sum)

## check PGG values for cancer pom replaced model-----------------
location.of.generated.files <- "~/Downloads/prostate-active-surveillance-vDaan/generated-files-cancerpom-replaced"
load(paste(location.of.generated.files,"IOP-data-shaping-work-space.RData", sep="/"))
(n <- dim(pt.data)[1])
#get observed PGGs
eta.data<-pt.data$true.pgg
obs.eta<-eta.data[!is.na(eta.data)]
obs.eta.mat<-cbind(as.numeric(obs.eta==1), as.numeric(obs.eta==2), as.numeric(obs.eta==3))
etahat<-read.csv(paste0(location.of.generated.files, "/jags-prediction-cancer_state-2022.csv"))
etahat<-as.matrix(etahat[,2:dim(etahat)[2]])
K<-4 
(B<-dim(etahat)[1]) 
(N<-length(obs.eta))
pred_eta <- matrix(nrow=(n-npat_cancer_known), ncol=K)
for(i in 1:(n-npat_cancer_known)){for(k in 1:K){
  pred_eta[i,k] <- sum(etahat[,i]==k)/B }}
alpha<-read.csv(paste0(location.of.generated.files, "/jags-prediction-pgg_intercept-2022.csv"))
alpha<-as.matrix(alpha[,2:dim(alpha)[2]])
alpha_med<-apply(alpha,2,median)
gamma<-read.csv(paste0(location.of.generated.files, "/jags-prediction-pgg_slope-2022.csv"))
gamma<-as.matrix(gamma[,2:dim(gamma)[2]])
gamma_med<-apply(gamma,2,median)
V.PGG.eta<-matrix(nrow=npat_pgg, ncol=(K-1))
for(j in 1:npat_pgg){
  if(pgg_patient_index_map[j] <= npat_cancer_known){
    for(k in 1:(K-1)){
      V.PGG.eta[j,k] <- as.numeric(eta.data[pgg_patient_index_map[j]]>k)}}
  else{
    for(k in 1:(K-1)){
      V.PGG.eta[j,k] <- 1-sum(pred_eta[(pgg_patient_index_map[j]-npat_cancer_known),1:k]) } } }
lin_pred1 <- alpha_med[1] - as.matrix(cbind(modmat_pgg, V.PGG.eta)) %*% as.vector(gamma_med)
lin_pred2 <- alpha_med[2] - as.matrix(cbind(modmat_pgg, V.PGG.eta)) %*% as.vector(gamma_med)
lin_pred3 <- alpha_med[3] - as.matrix(cbind(modmat_pgg, V.PGG.eta)) %*% as.vector(gamma_med)

exp1_canrep <- expit(lin_pred1)
exp2_canrep <- expit(lin_pred2) - exp1_canrep
exp3_canrep <- expit(lin_pred3) - exp2_canrep - exp1_canrep
exp4_canrep <- 1- exp1_canrep-exp2_canrep-exp3_canrep
exp_canrep <- cbind(exp1_canrep,exp2_canrep,exp3_canrep,exp4_canrep)
apply(exp_canrep, 2, sum)

## check PGG values for original-----------------
location.of.generated.files <- "~/Downloads/prostate-active-surveillance-vDaan/generated-files"
load(paste(location.of.generated.files,"IOP-data-shaping-work-space.RData", sep="/"))
source(paste(location.of.r.scripts,"data-prep-for-jags.R",sep="/"))
(n <- dim(pt.data)[1])
eta.data<-pt.data$true.pgg
obs.eta<-eta.data[!is.na(eta.data)]
obs.eta.mat<-cbind(as.numeric(obs.eta>1), as.numeric(obs.eta>2), as.numeric(obs.eta>3))
etahat<-read.csv(paste0(location.of.generated.files, "/jags-prediction-eta.hat-2021.csv"))
etahat<-as.matrix(etahat[,2:dim(etahat)[2]])
K<-4 #number of classes
(B<-dim(etahat)[1]) #length of saved posterior chain
(N<-length(obs.eta))
pred_eta <- matrix(nrow=(n-n_eta_known), ncol=K)
for(i in 1:(n-n_eta_known)){for(k in 1:K){
  pred_eta[i,k] <- sum(etahat[,i]==k)/B }}

alpha<-read.csv(paste0(location.of.generated.files, "/jags-prediction-alpha-2021.csv"))
alpha<-as.matrix(alpha[,2:dim(alpha)[2]])
alpha_med<-apply(alpha,2,median)
gamma<-read.csv(paste0(location.of.generated.files, "/jags-prediction-gamma.PGG-2021.csv"))
gamma<-as.matrix(gamma[,2:dim(gamma)[2]])
gamma_med<-apply(gamma,2,median)
V.PGG.eta<-matrix(nrow=n_pgg, ncol=(K-1))
for(j in 1:n_pgg){
  if(subj_pgg[j] <= n_eta_known){
    for(k in 1:(K-1)){
      V.PGG.eta[j,k] <- as.numeric(eta.data[subj_pgg[j]]>k)}}
  else{
    for(k in 1:(K-1)){
      V.PGG.eta[j,k] <- 1-sum(pred_eta[(subj_pgg[j]-n_eta_known),1:k]) } } }
lin_pred1 <- alpha_med[1] - as.matrix(cbind(V.PGG.data, V.PGG.eta)) %*% as.vector(gamma_med)
lin_pred2 <- alpha_med[2] - as.matrix(cbind(V.PGG.data, V.PGG.eta)) %*% as.vector(gamma_med)
lin_pred3 <- alpha_med[3] - as.matrix(cbind(V.PGG.data, V.PGG.eta)) %*% as.vector(gamma_med)

exp1_o <- expit(lin_pred1)
exp2_o <- expit(lin_pred2) - exp1_o
exp3_o <- expit(lin_pred3) - exp2_o - exp1_o
exp4_o <- 1- exp1_o-exp2_o-exp3_o
exp_o <- cbind(exp1_o,exp2_o,exp3_o,exp4_o)
apply(exp_o, 2, sum)

obs.pgg.mat<-cbind(as.numeric(pgg_data==1), as.numeric(pgg_data==2), 
                   as.numeric(pgg_data==3), as.numeric(pgg_data==4))
exp_o_tot <- apply(exp_o, 2, sum)
exp_canrep_tot <- apply(exp_canrep, 2, sum)
exp_seq_tot <- apply(exp_seq, 2, sum)
obs_tot <- apply(obs.pgg.mat, 2, sum)

get_chisq <- function(obs, exp){
  (obs[1]-exp[1])^2/exp[1] + (obs[2]-exp[2])^2/exp[2] + 
    (obs[3]-exp[3])^2/exp[3] + (obs[4]-exp[4])^2/exp[4]
}
get_chisq(obs_tot, exp_o_tot)
get_chisq(obs_tot, exp_canrep_tot)
get_chisq(obs_tot, exp_seq_tot)

### check AUC for latent cancer state (eta for patients with surgery) ------
#### Yate's model
rho_int <- read_csv("generated-files/jags-prediction-rho_int-2021.csv")
rho_coef <- read_csv("generated-files/jags-prediction-rho_coef-2021.csv")
rho_int <- rho_int[,-1]
rho_coef <- rho_coef[,-1]
rho_int_mean <-apply(rho_int, 2, mean) 
rho_coef_mean <-matrix(apply(rho_coef, 2, mean))

source(paste("R-scripts","data-prep-for-jags.R",sep="/"))
id_eta_obs <- pt.data$clinical_PTnum[!is.na(pt.data$true.pgg)]
eta_obs <- pt.data$true.pgg[!is.na(pt.data$true.pgg)]
inx_eta_obs <- which(!is.na(pt.data$true.pgg))
V.ETA.data.obs <- V.ETA.data[inx_eta_obs,]

linpred1 <- rho_int_mean[1] -  V.ETA.data.obs %*% rho_coef_mean
p1_eta <- exp(linpred1)/(1+exp(linpred1))

linpred2 <- rho_int_mean[2] -  V.ETA.data.obs %*% rho_coef_mean
cum_p2_eta <- exp(linpred2)/(1+exp(linpred2))
p2_eta <- cum_p2_eta - p1_eta

linpred3 <- rho_int_mean[3] -  V.ETA.data.obs %*% rho_coef_mean
cum_p3_eta <- exp(linpred3)/(1+exp(linpred3))
p3_eta <- cum_p3_eta - p1_eta - p2_eta

p4_eta <- 1-cum_p3_eta

obs1_eta <- as.numeric(eta_obs <= 2)
auc(obs1_eta, as.numeric(p1_eta+p2_eta))

obs1 <- sum(eta_obs==1); exp1 <- sum(p1_eta)
obs2 <- sum(eta_obs==2); exp2 <- sum(p2_eta)
obs3 <- sum(eta_obs==3); exp3 <- sum(p3_eta)
obs4 <- sum(eta_obs==4); exp4 <- sum(p4_eta)
obs_eta_original <- c(obs1, obs2, obs3, obs4)
exp_eta_original <- c(exp1, exp2, exp3, exp4)
get_chisq(obs_eta_original, exp_eta_original)

#### sequential replaced model
cancer_int1 <- read_csv("generated-files-seq/jags-prediction-cancer_int1-2022.csv")
cancer_int2 <- read_csv("generated-files-seq/jags-prediction-cancer_int2-2022.csv")
cancer_int3 <- read_csv("generated-files-seq/jags-prediction-cancer_int3-2022.csv")

cancer_slope1 <- read_csv("generated-files-seq/jags-prediction-cancer_slope1-2022.csv")
cancer_slope2 <- read_csv("generated-files-seq/jags-prediction-cancer_slope2-2022.csv")
cancer_slope3 <- read_csv("generated-files-seq/jags-prediction-cancer_slope3-2022.csv")

cancer_int1 <- cancer_int1[,-1];cancer_int2 <- cancer_int2[,-1];cancer_int3 <- cancer_int3[,-1]
cancer_slope1 <- cancer_slope1[,-1];cancer_slope2 <- cancer_slope2[,-1];cancer_slope3 <- cancer_slope3[,-1]

int1_mean <-apply(cancer_int1, 2, mean) 
int2_mean <-apply(cancer_int2, 2, mean) 
int3_mean <-apply(cancer_int3, 2, mean) 

slope1_mean <-apply(cancer_slope1, 2, mean) 
slope2_mean <-apply(cancer_slope2, 2, mean) 
slope3_mean <-apply(cancer_slope3, 2, mean) 

source(paste("R-scripts","data-prep-for-jags-reform.R",sep="/"))


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
get_chisq(obs_eta_original, exp_eta_seq)

obs1_eta <- as.numeric(eta_obs > 1)
auc(obs1_eta, 1-p1_eta_seq)


