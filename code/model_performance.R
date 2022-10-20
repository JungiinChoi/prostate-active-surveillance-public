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
alpha<-read.csv(paste0(location.of.generated.files, "/jags-prediction-pgg_intercept-2022.csv"))
alpha<-as.matrix(alpha[,2:dim(alpha)[2]])
# posterior estimate
alpha_med<-apply(alpha,2,median)

#gamma<-read.csv(paste0(location.of.generated.files, "/jags-prediction-gamma.PGG-2021.csv"))
gamma<-read.csv(paste0(location.of.generated.files, "/jags-prediction-pgg_slope-2022.csv"))
gamma<-as.matrix(gamma[,2:dim(gamma)[2]])
# posterior estimate
gamma_med<-apply(gamma,2,median)

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
#linear predictor
lin_pred <- alpha_med[1] - as.matrix(cbind(modmat_pgg, V.PGG.eta)) %*% as.vector(gamma_med)
#on probability scale
p_rc <- 1 - expit(lin_pred)

#observed outcome
RC <- as.numeric(pgg_data>1)

#loess fit for calibration plot with confidence intervals
loe.rc<-loess(RC~p_rc)
pred.loe<-predict(loe.rc, seq(0, 1, 0.01), se=T)
pred.loe.low <- pred.loe$fit + qnorm(0.025)*pred.loe$se.fit
pred.loe.hi <- pred.loe$fit + qnorm(0.975)*pred.loe$se.fit

#AUC for RC predictions
my.auc<-performance(prediction(p_rc, RC), "auc")@y.values[[1]]  #0.887
my.ci<-ci.auc(response=RC, predictor=p_rc, method="bootstrap") #takes a long time! #0.85, 0.89


#ROC
plot(c(0,1),c(0,1), type="n", xlab="False Positive Rate", 
     ylab="True Positive Rate")#, xaxt="n", yaxt="n")
roc.rc<-performance(prediction(p_rc, RC),"tpr","fpr")
lines(roc.rc@y.values[[1]]~roc.rc@x.values[[1]], lwd=3, col="blue") 
legend("bottomright", bty="n",
       legend=paste0("AUC = ", round(my.auc,2), " (95% CI = ", 
                     round(my.ci[1], 2), ", ",
                     round(my.ci[2], 2), ")"))


## for PGG >2 ------------------
#linear predictor
lin_pred2 <- alpha_med[2] - as.matrix(cbind(modmat_pgg, V.PGG.eta)) %*% as.vector(gamma_med)
#on probability scale
p_rc2 <- 1 - expit(lin_pred2)

#observed outcome
RC2 <- as.numeric(pgg_data>2)

#AUC for RC predictions
my.auc2<-performance(prediction(p_rc2, RC2), "auc")@y.values[[1]]  #0.887
my.ci2<-ci.auc(response=RC2, predictor=p_rc2, method="bootstrap") #takes a long time! #0.85, 0.89

## for PGG >3 ------------------
#linear predictor
lin_pred3 <- alpha_med[3] - as.matrix(cbind(modmat_pgg, V.PGG.eta)) %*% as.vector(gamma_med)
#on probability scale
p_rc3 <- 1 - expit(lin_pred3)

#observed outcome
RC3 <- as.numeric(pgg_data>3)

#AUC for RC predictions
my.auc2<-performance(prediction(p_rc3, RC3), "auc")@y.values[[1]]  #0.887
my.ci3<-ci.auc(response=RC3, predictor=p_rc3, method="bootstrap") #takes a long time! #0.85, 0.89
