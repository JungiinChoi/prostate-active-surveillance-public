# Yates Coley
# rycoley@gmail.com
# 2019 January 08
# jags model definition, updated for 10-fold CV

cat("model {


###PRIORS FOR LATENT CLASS MODEL

#flat priors for proportional odds model intercepts and coefficients with ordinal eta as outcome
#thresholds/intercept 
for(k in 1:(K-1)){rho_int0[k] ~ dnorm(0,0.01)}
rho_int[1:(K-1)] <- sort(rho_int0)
#proportional odds coefficients
for(index in 1:d.V.ETA){rho_coef[index] ~ dnorm(0,0.01)}


###PRIORS FOR PSA MIXED MODEL

#model correlated random effects distribution
for (index in 1:d.Z) {
	xi[index]~dunif(0,100)
	for(k in 1:K.bin){
		mu_raw[index,k]~dnorm(0, 0.01) 
		mu[index,k]<-xi[index] *mu_raw[index,k]}  }
for(k in 1:K.bin){
	mu_int[k] <- mu[1,k] 
	mu_slope[k] <- mu[2,k]}

#same covariance matrix (Sigma_B) across latent classes
Tau_B_raw ~ dwish(I_d.Z[,], (d.Z+1)) 
Sigma_B_raw[1:d.Z, 1:d.Z] <- inverse(Tau_B_raw[1:d.Z, 1:d.Z])	
for (index in 1:d.Z){
		sigma[index]<-xi[index]*sqrt(Sigma_B_raw[index,index]) }
sigma_int <- sigma[1] 
sigma_slope <- sigma[2] 
rho_int_slope <- Sigma_B_raw[1,2]/sqrt(Sigma_B_raw[1,1] * Sigma_B_raw[2,2])
cov_int_slope <- rho_int_slope*sigma_int*sigma_slope 

#residual variance, independent of correlated random effects, same across classes
sigma_res ~ dunif(0, 1) 
tau_res <- pow(sigma_res,-2)

#fixed effects
for(index in 1:d.X){
	beta[index] ~ dnorm(0,0.01)}



###PRIORS FOR BIOPSY OUTCOME MODEL
#additional element(s) in parameter vector are coefficients for class membership eta

##proportional odds regression for biopsy grade
#thresholds/intercept 
for(k in 1:(K-1)){alpha0[k] ~ dnorm(0,0.01)}
alpha[1:(K-1)] <- sort(alpha0)
#proportional odds coefficients
for(index in 1:d.V.PGG){gamma.PGG[index] ~ dnorm(0,0.01)}
for(index in (d.V.PGG+1):(d.V.PGG+(K-1))){gamma.PGG[index] ~ dnorm(0,0.01)T(-1,)}


###LIKELIHOOD

##latent variable for true cancer state
##proportional odds regression 
for(i in 1:n){
  lin.pred.eta[i] <- inprod(rho_coef[1:d.V.ETA], V.ETA[i,1:d.V.ETA])
  logit(cuml.p_eta[i,1]) <- rho_int[1] - lin.pred.eta[i]
  p_eta[i,1] <- cuml.p_eta[i,1]
  for(k in 2:(K-1)){
    logit(cuml.p_eta[i,k]) <- rho_int[k] - lin.pred.eta[i]
    p_eta[i,k] <- cuml.p_eta[i,k] - cuml.p_eta[i,(k-1)]}
  p_eta[i,K] <- 1 - cuml.p_eta[i,(K-1)] }
    
for(i in 1:n_eta_known){
  eta.data[i] ~ dcat(p_eta[i,1:K])
  eta[i] <- eta.data[i]
  eta.bin[i] <- step(eta[i]-2) + 1} #this is for those with path reports from SURG, eta known 
    
for(i in (n_eta_known+1):n){
  eta.hat[(i-n_eta_known)] ~ dcat(p_eta[i,1:K])
  eta[i] <- eta.hat[(i-n_eta_known)]
  eta.bin[i] <- step(eta[i]-2) + 1}  #for those without SURG
    
eta.track[1:n_mask] <- eta.hat[1:n_mask]

##linear mixed effects model for PSA 
#generate random intercept and slope for individual given latent class
for (i in 1:n) {
	B_raw[i,1:d.Z] ~ dmnorm(mu_raw[1:d.Z, eta.bin[i]], Tau_B_raw[1:d.Z, 1:d.Z])
	for(index in 1:d.Z){b.vec[i,index] <- xi[index]*B_raw[i,index]} }

#fit LME
for(j in 1:n_obs_psa){ 
	mu_obs_psa[j] <- inprod(b.vec[subj_psa[j],1:d.Z], Z[j,1:d.Z])  + inprod(beta[1:d.X], X[j,1:d.X]) 
	Y[j] ~ dnorm(mu_obs_psa[j], tau_res) }



### BIOPSY OUTCOMES

##proportional odds regression for biopsy upgrading 	
for(j in 1:n_pgg){
	lin.pred[j] <- inprod(gamma.PGG[1:d.V.PGG], V.PGG[j,1:d.V.PGG]) + gamma.PGG[(d.V.PGG+1)]*step(eta[subj_pgg[j]]-2) + gamma.PGG[(d.V.PGG+2)]*step(eta[subj_pgg[j]]-3) + gamma.PGG[(d.V.PGG+3)]*step(eta[subj_pgg[j]]-4)
	logit(cuml.p_rc[j,1]) <- alpha[1] - lin.pred[j]
	p_rc[j,1] <- cuml.p_rc[j,1]
	for(k in 2:(K-1)){
		logit(cuml.p_rc[j,k]) <- alpha[k] - lin.pred[j]	
		p_rc[j,k] <- cuml.p_rc[j,k] - cuml.p_rc[j,(k-1)] }
	p_rc[j,K] <- 1 - cuml.p_rc[j,(K-1)]
	PGG[j] ~ dcat(p_rc[j,1:K]) }

	
 }", fill=TRUE, file=paste(location.of.r.scripts,"cv-JAGS-prediction-model.txt", sep="/"))