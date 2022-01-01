# R.Y. Coley
# rycoley@gmail.com
# Updated 2019 January 07
# This code makes calibration plots for biopsy recalibration predictions and saves calibration objects


#### WORKFLOW
### 1. Define helper functions
### 2. Load model results
### 3. Get fitted biopsy reclaassification risk 
###
###



### 1. Define helper functions

# define inverse logit function
expit <- function(x)
{return(exp(x)/(1+exp(x)))}


### 2. Load model results

## get eta predictions
etahat<-read.csv(paste0(location.of.generated.files, "/jags-prediction-eta.hat-1.csv"))
etahat<-as.matrix(etahat[,2:dim(etahat)[2]])
for(i in 2:5){
  res<-read.csv(paste0(location.of.generated.files, "/jags-prediction-eta.hat-",i,".csv"))
  etahat<-rbind(etahat,res[,2:dim(res)[2]])}

#mean eta predictions
pred_eta <- matrix(nrow=(n-n_eta_known), ncol=K)
for(i in 1:(n-n_eta_known)){for(k in 1:K){
  pred_eta[i,k] <- sum(etahat[,i]==k)/B }}


#length of saved posterior chain
(B<-dim(etahat)[1])

## get parameter estimates in biopsy outcome proportional odds model
alpha<-read.csv(paste0(location.of.generated.files, "/jags-prediction-alpha-1.csv"))
alpha<-as.matrix(alpha[,2:dim(alpha)[2]])
for(i in 2:5){
  res<-read.csv(paste0(location.of.generated.files, "/jags-prediction-alpha-",i,".csv"))
  alpha<-rbind(alpha,res[,2:dim(res)[2]])}
# posterior estimate
alpha_med<-apply(alpha,2,median)

gamma<-read.csv(paste0(location.of.generated.files, "/jags-prediction-gamma.PGG-1.csv"))
gamma<-as.matrix(gamma[,2:dim(gamma)[2]])
for(i in 2:5){
  res<-read.csv(paste0(location.of.generated.files, "/jags-prediction-gamma.PGG-",i,".csv"))
  gamma<-rbind(gamma,res[,2:dim(res)[2]])}
# posterior estimate
gamma_med<-apply(gamma,2,median)




### 3. Get fitted biopsy reclaassification risk 

## update covariate matrix for biopsy model with eta predictions
V.PGG.eta<-matrix(nrow=n_pgg, ncol=(K-1))
for(j in 1:n_pgg){
  if(subj_pgg[j] <= n_eta_known){
    for(k in 1:(K-1)){
      V.PGG.eta[j,k] <- as.numeric(eta.data[subj_pgg[j]]>k)}}else{
        for(k in 1:(K-1)){
          V.PGG.eta[j,k] <- 1-sum(pred_eta[(subj_pgg[j]-n_eta_known),1:k]) } } }


## get estimated probability of upgrading

#linear predictor
lin_pred <- alpha_med[1] - as.matrix(cbind(V.PGG.data, V.PGG.eta)) %*% as.vector(gamma_med)
#on probability scale
p_rc <- 1 - expit(lin_pred)

#observed outcome
RC <- as.numeric(PGG>1)

#loess fit for calibration plot with confidence intervals
loe.rc<-loess(RC~p_rc)
pred.loe<-predict(loe.rc, seq(0, 1, 0.01), se=T)
pred.loe.low <- pred.loe$fit + qnorm(0.025)*pred.loe$se.fit
pred.loe.hi <- pred.loe$fit + qnorm(0.975)*pred.loe$se.fit

#what is the AUC for RC predictions?
my.auc<-performance(prediction(p_rc, RC), "auc")@y.values[[1]]  #0.887
my.ci<-ci.auc(response=RC, predictor=p_rc, method="bootstrap") #takes a long time! #0.85, 0.89


#print calibration curve and ROC curve on same pdf
pdf(paste0(location.of.generated.files,"/predictive-accuracy-rc.pdf"), 
    width=10, height=5)

par(mfrow=c(1,2))

#calibration plot
plot(RC~p_rc, cex=0.25, xlab="Predicted risk of upgrading",
     ylab="Observed proportion of biopsies with upgrading")
lines(c(0,1), c(0,1), col="black", lwd=3, lty="dotted")

lines(pred.loe$fit~seq(0, 1, 0.01), col="navy", lwd=3)
lines(pred.loe.low~seq(0, 1, 0.01), col="navy", lwd=3, lty="dashed")
lines(pred.loe.hi~seq(0, 1, 0.01), col="navy", lwd=3, lty="dashed")

#ROC
plot(c(0,1),c(0,1), type="n", xlab="False Positive Rate", 
     ylab="True Positive Rate")#, xaxt="n", yaxt="n")
roc.rc<-performance(prediction(p_rc, RC),"tpr","fpr")
lines(roc.rc@y.values[[1]]~roc.rc@x.values[[1]], lwd=3, col="blue") 
legend("bottomright", bty="n",
       legend=paste0("AUC = ", round(my.auc,2), " (95% CI = ", 
                     round(my.ci[1], 2), ", ",
                     round(my.ci[2], 2), ")"))
dev.off()