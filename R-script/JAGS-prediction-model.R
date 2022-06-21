#UPDATED 10/18/17
#Annotations updated 12/21/2017
#### IOP components and additinal biopsy outcomes (NPC, MPC, LAT) are commented out because they do not improve model estimation

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



###PRIORS FOR BIOPSY OUTCOMES AND SURGERY MODEL
#additional element(s) in each parameter vector is coefficient for class membership eta

##biopsy received
#for(index in 1:(d.U.BX+1)){nu.BX[index] ~ dnorm(0,0.01)}

#negative binomial regression for NPC
#regression coefficients
#for(index in 1:d.V.NPC){beta.NPC[index] ~ dnorm(0, 0.01)}
#beta.NPC[(d.V.NPC+1)] ~ dnorm(0, 0.01)T(-1,)
#variance
#r.NPC ~ dgamma(0.01, 0.01)

##proportional odds regression for biopsy grade
#thresholds/intercept
for(k in 1:(K-1)){alpha0[k] ~ dnorm(0,0.01)}
alpha[1:(K-1)] <- sort(alpha0)
#proportional odds coefficients
for(index in 1:d.V.PGG){gamma.PGG[index] ~ dnorm(0,0.01)}
#for(index in (d.V.PGG+1):(d.V.PGG+(K-1))){gamma.PGG[index] ~ dnorm(0,0.01)T(-1,)}
gamma.PGG[d.V.PGG+1] ~ dnorm(3.25,0.01)T(-1,) # truncate at left at -1
gamma.PGG[d.V.PGG+2] ~ dnorm(4.58,0.01)T(-1,)
gamma.PGG[d.V.PGG+3] ~ dnorm(5.84,0.01)T(-1,)

#negative binomial regression for MPC
#regression coefficients
#for(index in 1:d.V.MPC){beta.MPC[index] ~ dnorm(0, 0.01)}
#beta.MPC[(d.V.MPC+1)] ~ dnorm(0, 0.01)T(-1,)
#variance
#r.MPC ~ dgamma(0.01, 0.01)

## Laterality logistic regression
#for(index in 1:d.V.LAT){beta.LAT[index] ~ dnorm(0,0.01)}
#beta.LAT[(d.V.LAT+1)] ~ dnorm(0, 0.01)T(-1,)

##surgery received
#includes interaction with prior grade RC
#for(index in 1:(d.W.SURG+2)){omega.SURG[index] ~ dnorm(0,0.01)}



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


##linear mixed effects model for PSA
#generate random intercept and slope for individual given latent class
for (i in 1:n) {
	B_raw[i,1:d.Z] ~ dmnorm(mu_raw[1:d.Z, eta.bin[i]], Tau_B_raw[1:d.Z, 1:d.Z])
	for(index in 1:d.Z){b.vec[i,index] <- xi[index]*B_raw[i,index]} }

#fit LME
for(j in 1:n_obs_psa){
	mu_obs_psa[j] <- inprod(b.vec[subj_psa[j],1:d.Z], Z[j,1:d.Z])  + inprod(beta[1:d.X], X[j,1:d.X])
	Y[j] ~ dnorm(mu_obs_psa[j], tau_res) }



### BIOPSY OUTCOMES AND SURGERY RECEIVED

##logistic regression for biopsy
#for(j in 1:n_bx){
#	logit(p_bx[j]) <-inprod(nu.BX[1:d.U.BX], U.BX[j,1:d.U.BX]) + nu.BX[(d.U.BX+1)]*equals(eta.bin[subj_bx[j]],2)
#	BX[j] ~ dbern(p_bx[j]) }


##negative binomial regression for NPC
#for(j in 1:n_npc){
##log(NCS.offset[j]) + removed in favor of ns rep
#  log(mu_npc[j]) <-  inprod(beta.NPC[1:d.V.NPC], V.NPC[j,1:d.V.NPC]) + beta.NPC[(d.V.NPC+1)]*step(eta[subj_npc[j]]-2)
### + beta.NPC[(d.V.NPC+2)]*step(eta[subj_npc[j]]-3) + beta.NPC[(d.V.NPC+3)]*step(eta[subj_npc[j]]-4)
#  p_npc[j] <- r.NPC / (r.NPC + mu_npc[j])
#  NPC[j] ~ dnegbin(p_npc[j], r.NPC) }


##proportional odds regression for biopsy upgrading
for(j in 1:n_pgg){
	lin.pred[j] <- inprod(gamma.PGG[1:d.V.PGG], V.PGG[j,1:d.V.PGG]) + 
	gamma.PGG[(d.V.PGG+1)]*step(eta[subj_pgg[j]]-2) + 
	gamma.PGG[(d.V.PGG+2)]*step(eta[subj_pgg[j]]-3) + 
	gamma.PGG[(d.V.PGG+3)]*step(eta[subj_pgg[j]]-4)
	logit(cuml.p_rc[j,1]) <- alpha[1] - lin.pred[j]
	p_rc[j,1] <- cuml.p_rc[j,1]
	for(k in 2:(K-1)){
		logit(cuml.p_rc[j,k]) <- alpha[k] - lin.pred[j]
		p_rc[j,k] <- cuml.p_rc[j,k] - cuml.p_rc[j,(k-1)] }
	p_rc[j,K] <- 1 - cuml.p_rc[j,(K-1)]
	PGG[j] ~ dcat(p_rc[j,1:K]) }


##negative binomial regression for MPC
#for(j in 1:n_mpc){
#  log(mu_mpc[j]) <- log(10) + inprod(beta.MPC[1:d.V.MPC], V.MPC[j,1:d.V.MPC]) + beta.MPC[(d.V.MPC+1)]*step(eta[subj_mpc[j]]-2)
###+ beta.MPC[(d.V.MPC+2)]*step(eta[subj_mpc[j]]-3) + beta.MPC[(d.V.MPC+3)]*step(eta[subj_mpc[j]]-4)
#  p_mpc[j] <- r.MPC / (r.MPC + mu_mpc[j])
#  MPC[j] ~ dnegbin(p_mpc[j], r.MPC) }


##logistic regression for bilateral cancer
#for(j in 1:n_lat){
#	logit(p_lat[j]) <-inprod(beta.LAT[1:d.V.LAT], V.LAT[j,1:d.V.LAT]) + beta.LAT[(d.V.LAT+1)]*step(eta[subj_lat[j]]-2)
###+ beta.LAT[(d.V.LAT+2)]*step(eta[subj_lat[j]]-3) + beta.LAT[(d.V.LAT+3)]*step(eta[subj_lat[j]]-4)
#	LAT[j] ~ dbern(p_lat[j]) }


##logistic regression for surgery
#for(j in 1:n_surg){
#	logit(p_surg[j]) <-inprod(omega.SURG[1:d.W.SURG], W.SURG[j,1:d.W.SURG]) + omega.SURG[(d.W.SURG+1)]*equals(eta.bin[subj_surg[j]],2)  + omega.SURG[(d.W.SURG+2)]*equals(eta.bin[subj_surg[j]],2)*W.SURG[j,(d.W.SURG-1)]
#	SURG[j] ~ dbern(p_surg[j]) }


 }", fill=TRUE, file=paste(location.of.r.scripts,"JAGS-prediction-model.txt", sep="/"))
