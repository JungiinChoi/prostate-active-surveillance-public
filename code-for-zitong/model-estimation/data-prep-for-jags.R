#Yates Coley
#rycoley@gmail.com
#2018 February 18, updated 2018 March 16, updated March 19
#This script will take tidied/shaped AS data and perform additional manipulation for the JAGS model
#I've put a lot of data checks in here as fail-safes, but there shouldn't be any missing values after running the previous script.


#### WORKFLOW
### 1. Format pt-level data for JAGS run 
### 2. Format PSA data for JAGS run 
### 3. Format biopsy result (reclassification) data for JAGS run 



### 0 .Load libraries
library("lme4")
library("splines")


### 1. Format pt-level data for JAGS run 
#number of patients
(n <- dim(pt.data)[1])

#list of known true states
eta.data <- pt.data$true.pgg[!is.na(pt.data$true.pgg)]

#number with known true state
(n_eta_known <- length(eta.data))


#covariate matrix for predictors of PGG
#### NEED TO ADD MAT FOR RANDOM EFFECTS MODEL
V.ETA.data <- as.matrix(cbind(pt.data$dx.age.std))
d.V.ETA<-dim(V.ETA.data)[2]



### 2. Format PSA data for JAGS run 
#number of observations
(n_obs_psa <- dim(psa.data)[1])

#observed PSA
Y <- psa.data$log.psa
summary(Y)
data.check(condition=as.logical(sum(is.na(Y))==0), message="Missing PSA values. Email Yates; she will check orginal script.")

#list of patients
subj_psa <- psa.data$subj
sum(is.na(subj_psa))
length(unique(subj_psa))

#covariate matrix for random effects
psa.data$years <- psa.data$days/365.25
Z.data <- as.matrix(cbind(rep(1,n_obs_psa), psa.data$years))
d.Z <- dim(Z.data)[2]
data.check(condition=as.logical(sum(is.na(Z.data))==0), message="Missing time since dx in PSA data. Email Yates; she will check orginal script.")
summary(psa.data$days)

#covariate matrix for fixed effects
psa.data$dx.age.std<-scale(psa.data$dx.age.std)
X.data <- as.matrix(cbind(psa.data$vol.std, psa.data$dx.age.std))
d.X <- dim(X.data)[2]
data.check(condition=as.logical(sum(is.na(X.data))==0), message="Missing volumes in PSA data. Email Yates; she will check orginal script.")
summary(psa.data$vol.std)
summary(psa.data$dx.age.std)

#lmer fit to get starting value for covariance parameter
mod_lmer<-lmer(log.psa~ vol.std + dx.age.std + (1+ years |subj), data=psa.data)
var_vec <- apply(coef(mod_lmer)$subj, 2, var)[1:d.Z]
names(var_vec)
#make sure order of variance estimates is correct
index.intercept<-c(1:2)[names(var_vec)=="(Intercept)"]
index.time<-c(1:2)[names(var_vec)=="years"]
var_vec <- c(var_vec[index.intercept], var_vec[index.time])



### 3. Format biopsy result (reclassification) data for JAGS run 
#number of biopsies with GS outcome
(n_pgg <- dim(bx.data)[1])

#GS outcomes for biopsies
PGG <- as.numeric(bx.data$bx.pgg)
table(PGG)
data.check(condition=as.logical(sum(is.na(PGG))==0), message="Missing biopsy results data. Email Yates; she will check orginal script.")

#unique subj identifier
subj_pgg <- bx.data$subj

#covariate matrix
bx.yr.knots <- attr(ns(bx.data$yr, 3), "knots")
bx.yr.bknots <- attr(ns(bx.data$yr, 3), "Boundary.knots")
bx.ncs.knots <- attr(ns(bx.data$num_cores_sampled, 2), "knots")
bx.ncs.bknots <- attr(ns(bx.data$num_cores_sampled, 2), "Boundary.knots")
bx.vol.knots <- attr(ns(bx.data$vol.std, 2), "knots")
bx.vol.bknots <- attr(ns(bx.data$vol.std, 2), "Boundary.knots")

V.PGG.data <- as.matrix(cbind(ns(bx.data$yr, 3), 
                              ns(bx.data$num_cores_sampled,2),
                              ns(bx.data$vol.std,2) ))
d.V.PGG <- dim(V.PGG.data)[2]
summary(V.PGG.data)
data.check(condition=as.logical(sum(is.na(V.PGG.data))==0), message="Missing biopsy data. Email Yates; she will check orginal script.")








