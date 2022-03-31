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

obs_tot <- data.frame(obs_tot);colnames(obs_tot) <- "N"; obs_tot$model <- "observed"
exp_o_tot <- data.frame(exp_o_tot); colnames(exp_o_tot) <- "N";exp_o_tot$model <- "exp_original"
exp_canrep_tot <- data.frame(exp_canrep_tot); colnames(exp_canrep_tot) <- "N";exp_canrep_tot$model <- "exp_cancer_replaced"
exp_seq_tot <- data.frame(exp_seq_tot); colnames(exp_seq_tot) <- "N";exp_seq_tot$model <- "exp_both_replaced"

obs_tot$outcome <- exp_o_tot$outcome <- c("pgg=1", "pgg=2", "pgg=3", "pgg=4")
exp_canrep_tot$outcome <- exp_seq_tot$outcome <- c("pgg=1", "pgg=2", "pgg=3", "pgg=4")

bardat <- rbind(obs_tot,exp_o_tot,exp_canrep_tot,exp_seq_tot)
bardat$model <- factor(bardat$model, levels = c("observed", "exp_original", "exp_cancer_replaced", "exp_both_replaced"))
p1 <- ggplot(data = bardat %>% filter(outcome == "pgg=1"), aes(outcome, N, fill = model))+
  geom_bar(stat ='identity', position = 'dodge')

p2 <- ggplot(data = bardat %>% filter(outcome != "pgg=1"), aes(outcome, N, fill = model))+
  geom_bar(stat ='identity', position = 'dodge')
library(gridExtra)
grid.arrange(p1, p2, nrow= 1)
