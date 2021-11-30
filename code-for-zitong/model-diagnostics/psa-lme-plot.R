# Yates Coley
# rycoley@gmail.com
# 2019 January 07
# This script visualizes PSA posterior estimates

#### WORKFLOW
### 1. Load individual random effect estimates
### 2. Shape other parameter estimates, individual predictions (previously loaded)
### 3. Create visualization



### 1. Load individual random effect estimates
# get individual's random effects
b.vec<-read.csv(paste0(location.of.generated.files, "/jags-prediction-b.vec-1.csv"))
b.int<-as.matrix(b.vec[,2:(n+1)])
b.slope<-as.matrix(b.vec[,(n+2):(2*n+1)])

b.int<-apply(b.int,2,mean)
b.slope<-apply(b.slope,2,mean)


### 2. Shape other parameter estimates, individual predictions (previously loaded)
##paramters of PSA model
#mean for each group
mu_int<-apply(mu_int,2,mean)
mu_slope<-apply(mu_slope,2,mean)
mu0<-c(mu_int[1], mu_slope[1])
mu1<-c(mu_int[2], mu_slope[2])

#variance estimates
var_int<-mean(sigma_int)^2
var_slope<-mean(sigma_slope)^2
cov_int_slope<-mean(cov_int_slope)
#make covariance matrix
(Sigma<-matrix(c(var_int, cov_int_slope, cov_int_slope, var_slope), nrow=2, ncol=2))


#get individual's true state predictions (or observations)
eta.category<-vector(length=n)
eta.category[!is.na(pt.data$true.pgg) & pt.data$true.pgg==1] <- "Indolent"
eta.category[!is.na(pt.data$true.pgg) & pt.data$true.pgg>1] <- "Agressive"

eta.pred<-vector(length=(n-n_eta_known))
for(i in 1:(n-n_eta_known)){
	eta.pred[i]<-mean(etahat[,i]>1)
	if(eta.pred[i]<=0.25){eta.category[(n_eta_known+i)]<-"0-25%"}
	if(eta.pred[i]>0.25 & eta.pred[i]<=0.5){eta.category[(n_eta_known+i)]<-"25-50%"}
	if(eta.pred[i]>0.5 & eta.pred[i]<=0.75){eta.category[(n_eta_known+i)]<-"50-75%"}
	if(eta.pred[i]>0.75){eta.category[(n_eta_known+i)]<-"75-100%"}
	}
table(eta.category)


### 3. Create visualization

#make plot
pdf(paste0(location.of.assessment.summaries, "/psa-lme-plot.pdf"), width=9, height=7)

plot(c(min(b.slope)-0.1, max(b.slope)+0.1)~c(min(b.int), max(b.int)), type="n", xlab="Intercept", ylab="Slope" )

ellipse(mu0, Sigma, alpha=0.05, lwd=3, col="green3", lty="dotted")
ellipse(mu0, Sigma, alpha=0.5, lwd=3, col="green3", lty="dashed")
ellipse(mu1, Sigma, alpha=0.05, lwd=3, col="red", lty="dotted")
ellipse(mu1, Sigma, alpha=0.5, lwd=3, col="red", lty="dashed")



for(i in 1:sum(eta.category=="Indolent")){
	points(b.slope[i]~b.int[i], pch=19, col="green3", cex=0.4)}
for(i in (sum(eta.category=="Indolent")+1):n_eta_known){
	points(b.slope[i]~b.int[i], pch=19, col="red", cex=0.4)}

for(i in (n_eta_known+1):n){
	if(eta.pred[(i-n_eta_known)]<=quantile(eta.pred,0.25)){points(b.slope[i]~b.int[i], pch=21, col="green3", cex=0.4)}
	else{ if(eta.pred[(i-n_eta_known)]>quantile(eta.pred,0.25) & eta.pred[(i-n_eta_known)]<=median(eta.pred)){points(b.slope[i]~b.int[i], pch=21, col="greenyellow", cex=0.4)}
		else{
	if(eta.pred[(i-n_eta_known)]<=quantile(eta.pred,0.75) & eta.pred[(i-n_eta_known)]>median(eta.pred)){points(b.slope[i]~b.int[i], pch=21, col="orange", cex=0.4)}
		else{
	if(eta.pred[(i-n_eta_known)]>quantile(eta.pred,0.75)){points(b.slope[i]~b.int[i], pch=21, col="red", cex=0.4)}}
	}}}


legend("topleft", title="P(PGG>1)", legend=c("0-25%", "26-50%", "51-75%", "75-100%"), col=c("green3", "greenyellow", "orange", "red"), pch=rep(21,4), bty="n")
legend(0.2,1.75, title="State Observed", legend=c("PGG=1", "PGG>1"), col=c("green3","red"), pch=rep(19,2), bty="n")

legend("bottomleft", title="Credible Intervals", legend=c("50%", "95%"), lty=c("dashed","dotted"), lwd=rep(3,2), bty="n")
dev.off()