
# Yates Coley
# rycoley@gmail.com
# 2019 January 07
# This code generates diagnostic plots for model convergence

#### WORKFLOW
### 1. Define helper functions
### 2. Convergence plot for latent class predictions
### 3. For parameter estimates: load posterior results, make convergence plots


### 1. Define helper functions

##define inverse logit function
expit<-function(x){return(exp(x)/(1+exp(x)))}

##define functions for each plot type 
#each color represents one starting seed/separate MCMC chain
color.post<-c("blue","red","green","orange","navy")

#trace plots for all chains
plot.post.trace<-function(par,nc,n.row, n.col){
	num.ea<-dim(par)[1]/nc
	par(mfrow=c(n.row,n.col), mar=c(2,2,1,1))
	for(i in 1:dim(par)[2]){
		plot(c(min(par[,i]),max(par[,i]))~c(1,num.ea), type="n", xlab="", ylab=paste(i))
		for(j in 1:nc){
			lines(par[((j-1)*num.ea+1):(j*num.ea),i], col=color.post[j])}}}
		
#posterior density plots for all chains
plot.post.dens<-function(par, nc, n.row, n.col){
	num.ea<-dim(par)[1]/nc
	par(mfrow=c(n.row,n.col), mar=c(2,2,1,1))
	for(i in 1:dim(par)[2]){
		plot(c(0, max(density(par[,i])$y)*1.1) ~ c(min(par[,i]),max(par[,i])), type="n", xlab=paste(i), ylab="Density")
		for(j in 1:nc){
			lines(density(par[((j-1)*num.ea+1):(j*num.ea),i]), col=color.post[j])}}}

#cumulative quantile plot for one chain
plot.post.cuml<-function(par, nc, n.row, n.col){
	num.ea<-dim(par)[1]/nc
	par(mfrow=c(n.row,n.col), mar=c(2,2,1,1))
	for(i in 1:dim(par)[2]){
		cumuplot(par[1:num.ea,i], ylab="Quantile Estimate", xlab="Iteration", auto.layout=F)}}
#I had difficulty getting separate cumulative plot lines on an existing plot
#		for(j in 1:nc){	lines(cumuplot(par[((j-1)*num.ea+1):(j*num.ea),i]), col=color.post[j]) }} }



### 2. Convergence plot for latent class predictions
etahat<-read.csv(paste(location.of.generated.files, "/jags-prediction-eta.hat-1.csv",sep=""))
etahat<-as.matrix(etahat[,2:dim(etahat)[2]])
for(i in 2:5){
	res<-read.csv(paste(location.of.generated.files, "/jags-prediction-eta.hat-",i,".csv",sep=""))
	etahat<-rbind(etahat,res[,2:dim(res)[2]])}
	
#single mean per person
eta.mean<-apply(etahat,2,mean)

#get mean per person per chain
(a<-dim(etahat)[1]/5*c(0:4)+1)
(b<-dim(etahat)[1]/5*c(1:5))
eta.mean.mat<-as.matrix(cbind(apply(etahat[a[1]:b[1],],2,mean), 
                              apply(etahat[a[2]:b[2],],2,mean), 
                              apply(etahat[a[3]:b[3],],2,mean), 
                              apply(etahat[a[4]:b[4],],2,mean), 
                              apply(etahat[a[5]:b[5],],2,mean)))

pdf(paste0(location.of.assessment.summaries, "/post-density-eta.pdf"))
plot(density(eta.mean.mat[,1]),type="l",col="blue",xlim=c(1,4), ylim=c(0,max(density(eta.mean.mat[,1])$y) + 0.2), lwd=2, main="Probablity Lethal PCa")
lines(density(eta.mean.mat[,2]),col="red",lwd=2)
lines(density(eta.mean.mat[,3]),col="green",lwd=2)
lines(density(eta.mean.mat[,4]),col="orange",lwd=2)
lines(density(eta.mean.mat[,5]),col="navy",lwd=2)
dev.off()



### 3. For parameter estimates: load posterior results, make convergence plots

## RHO intercept for proportional odds model
rho_int<-read.csv(paste(location.of.generated.files, "/jags-prediction-rho_int-1.csv",sep=""))
rho_int<-as.matrix(rho_int[,2:dim(rho_int)[2]])
for(i in 2:5){
	res<-read.csv(paste(location.of.generated.files, "/jags-prediction-rho_int-",i,".csv",sep=""))
  rho_int<-rbind(rho_int,res[,2:dim(res)[2]])}

pdf(paste0(location.of.assessment.summaries, "/post-trace-rho_int.pdf"))
plot.post.trace(as.matrix(rho_int),5,2,2) 
dev.off()

pdf(paste0(location.of.assessment.summaries, "/post-dens-rho_int.pdf"))
plot.post.dens(as.matrix(rho_int),5,2,2)
dev.off()

pdf(paste0(location.of.assessment.summaries, "/post-cuml-rho_int.pdf"))
plot.post.cuml(as.matrix(rho_int),5,2,2)
dev.off()


## RHO coefficients for proportional odds model

rho_coef<-read.csv(paste(location.of.generated.files, "/jags-prediction-rho_coef-1.csv",sep=""))
rho_coef<-as.matrix(rho_coef[,2])
for(i in 2:5){
  res<-read.csv(paste(location.of.generated.files, "/jags-prediction-rho_coef-",i,".csv",sep=""))
  rho_coef<-c(rho_coef,res[,2])}

pdf(paste0(location.of.assessment.summaries, "/post-trace-rho_coef.pdf"))
plot.post.trace(as.matrix(rho_coef),5,2,2) 
dev.off()

pdf(paste0(location.of.assessment.summaries, "/post-dens-rho_coef.pdf"))
plot.post.dens(as.matrix(rho_coef),5,2,2)
dev.off()

pdf(paste0(location.of.assessment.summaries, "/post-cuml-rho_coef.pdf"))
plot.post.cuml(as.matrix(rho_coef),5,2,2)
dev.off()





### biopsy grade reclassification model
##intercept for proportional odds model
alpha<-read.csv(paste(location.of.generated.files, "/jags-prediction-alpha-1.csv",sep=""))
dim(alpha)
alpha<-as.matrix(alpha[,2:dim(alpha)[2]])
for(i in 2:5){
	res<-read.csv(paste(location.of.generated.files, "/jags-prediction-alpha-",i,".csv",sep=""))
	alpha<-rbind(alpha,res[,2:dim(res)[2]])}

pdf(paste0(location.of.assessment.summaries, "/post-trace-alpha.pdf"))
plot.post.trace(as.matrix(alpha),5,2,2)
dev.off()

pdf(paste0(location.of.assessment.summaries, "/post-dens-alpha.pdf"))
plot.post.dens(as.matrix(alpha),5,2,2)
dev.off()

pdf(paste0(location.of.assessment.summaries, "/post-cuml-alpha.pdf"))
plot.post.cuml(as.matrix(alpha),5,2,2) 
dev.off()

##coefficients for proportional odds model
gamma.PGG<-read.csv(paste(location.of.generated.files, "/jags-prediction-gamma.PGG-1.csv",sep=""))
dim(gamma.PGG)
gamma.PGG<-as.matrix(gamma.PGG[,2:dim(gamma.PGG)[2]])
for(i in 2:5){
	res<-read.csv(paste(location.of.generated.files, "/jags-prediction-gamma.PGG-",i,".csv",sep=""))
	gamma.PGG<-rbind(gamma.PGG,res[,2:dim(res)[2]])}

#apply(gamma.PGG,2,summary) 
d_pgg<-dim(gamma.PGG)[2]

pdf(paste0(location.of.assessment.summaries, "/post-trace-gamma_rc.pdf"))
plot.post.trace(as.matrix(gamma.PGG),5,(d_pgg/3),3)
dev.off()

pdf(paste0(location.of.assessment.summaries, "/post-dens-gamma_rc.pdf"))
plot.post.dens(as.matrix(gamma.PGG),5,(d_pgg/3),3)
dev.off()
#here, we see that the trace plots for each chain show convergence and that results are similar across sampling chains

pdf(paste0(location.of.assessment.summaries, "/post-cuml-gamma_rc.pdf"))
plot.post.cuml(as.matrix(gamma.PGG),5,(d_pgg/3),3) 
dev.off()




###PSA model

##fixed effect for age at diagnosis and prostate volume
beta<-read.csv(paste(location.of.generated.files, "/jags-prediction-beta-1.csv",sep=""))
dim(beta)
beta<-as.matrix(beta[,2:3])
for(i in 2:5){
	res<-read.csv(paste(location.of.generated.files, "/jags-prediction-beta-",i,".csv",sep=""))
	beta<-rbind(beta, res[,2:3])}
apply(beta,2,summary)

pdf(paste0(location.of.assessment.summaries, "/post-trace-beta.pdf"))
plot.post.trace(as.matrix(beta),5,2,2)
dev.off()

pdf(paste0(location.of.assessment.summaries, "/post-dens-beta.pdf"))
plot.post.dens(as.matrix(beta),5,2,2)
dev.off()

pdf(paste0(location.of.assessment.summaries, "/post-cuml-beta.pdf"))
plot.post.cuml(as.matrix(beta),5,2,2) 
dev.off()


#mean random intercept for each class
mu_int<-read.csv(paste(location.of.generated.files, "/jags-prediction-mu_int-1.csv",sep=""))
mu_int<-as.matrix(mu_int[,2:3])
for(i in 2:5){
	res<-read.csv(paste(location.of.generated.files, "/jags-prediction-mu_int-",i,".csv",sep=""))
	mu_int<-rbind(mu_int,res[,2:3])}
apply(mu_int,2,summary)

#mean random slope for each class
mu_slope<-read.csv(paste(location.of.generated.files, "/jags-prediction-mu_slope-1.csv",sep=""))
mu_slope<-as.matrix(mu_slope[,2:3])
for(i in 2:5){
	res<-read.csv(paste(location.of.generated.files, "/jags-prediction-mu_slope-",i,".csv",sep=""))
	mu_slope<-rbind(mu_slope,res[,2:3])}
apply(mu_slope,2,summary)

pdf(paste0(location.of.assessment.summaries, "/post-trace-mu.pdf"))
plot.post.trace(cbind(mu_int,mu_slope),5,2,2)
dev.off()

pdf(paste0(location.of.assessment.summaries, "/post-dens-mu.pdf"))
plot.post.dens(cbind(mu_int,mu_slope),5,2,2)
dev.off()
#here, we see that the trace plots for each chain show convergence and that results are similar across sampling chains

pdf(paste0(location.of.assessment.summaries, "/post-cuml-mu.pdf"))
plot.post.cuml(cbind(mu_int, mu_slope),5,2,2) 
dev.off()


#standard deviation for random intercepts
sigma_int<-read.csv(paste(location.of.generated.files, "/jags-prediction-sigma_int-1.csv",sep=""))
sigma_int<-as.matrix(sigma_int[,2])
for(i in 2:5){
	res<-read.csv(paste(location.of.generated.files, "/jags-prediction-sigma_int-",i,".csv",sep=""))
	sigma_int<-c(sigma_int,res[,2])}
summary(sigma_int)


#sandard deviation for random slope
sigma_slope<-read.csv(paste(location.of.generated.files, "/jags-prediction-sigma_slope-1.csv",sep=""))
sigma_slope<-as.matrix(sigma_slope[,2])
for(i in 2:5){
	res<-read.csv(paste(location.of.generated.files, "/jags-prediction-sigma_slope-",i,".csv",sep=""))
	sigma_slope<-c(sigma_slope,res[,2])}
summary(sigma_slope)

#covariance of random effects
cov_int_slope<-read.csv(paste(location.of.generated.files, "/jags-prediction-cov_int_slope-1.csv",sep=""))
cov_int_slope<-as.matrix(cov_int_slope[,2])
for(i in 2:5){
	res<-read.csv(paste(location.of.generated.files, "/jags-prediction-cov_int_slope-",i,".csv",sep=""))
	cov_int_slope<-c(cov_int_slope,res[,2])}
summary(cov_int_slope)

#residual variance for random effects model
sigma_res<-read.csv(paste(location.of.generated.files, "/jags-prediction-sigma_res-1.csv",sep=""))
sigma_res<-as.matrix(sigma_res[,2])
for(i in 2:5){
	res<-read.csv(paste(location.of.generated.files, "/jags-prediction-sigma_res-",i,".csv",sep=""))
	sigma_res<-c(sigma_res,res[,2])}
summary(sigma_res)



pdf(paste0(location.of.assessment.summaries, "/post-trace-cov-and-sigma.pdf"))
plot.post.trace(cbind(sigma_int, sigma_slope, cov_int_slope, sigma_res),5,2,2)
dev.off()

pdf(paste0(location.of.assessment.summaries, "/post-dens-cov-and-sigma.pdf"))
plot.post.dens(cbind(sigma_int, sigma_slope, cov_int_slope, sigma_res),5,2,2)
dev.off()


pdf(paste0(location.of.assessment.summaries, "/post-cuml-cov-and-sigma.pdf"))
plot.post.cuml(cbind(sigma_int, sigma_slope, cov_int_slope, sigma_res),5,2,2) 
dev.off()




