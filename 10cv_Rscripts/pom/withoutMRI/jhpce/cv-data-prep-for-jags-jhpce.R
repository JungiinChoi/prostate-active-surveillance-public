########### SEED/GROUP TO MASK IS EITHER SET HERE OR IN EXAMPLE-CV-TO-RUN.R
#### This should be adjusted by the user depending on how seed, masking set 
# SEED <- as.numeric(Sys.getenv("SGE_TASK_ID"))
# # For GAP3 analysis, make this seed=1,...,10 and define grouping variable with values 1,...,10 
# to.mask<- SEED
###########



# Yates Coley
# rycoley@gmail.com
# 2019 January 08, modified 2019 January 30
# This script will take tidied/shaped AS data and perform additional manipulation for the JAGS model
# I've put a lot of data checks in here as fail-safes, but there shouldn't be any missing values after running the previous script.
# This updated from data-prep-for-jags.R to accomodated 10-fold CV 

#### WORKFLOW
### 1. Format pt-level data for JAGS run
#      Include removing true state observation for CV leave-one-out fold
### 2. Format PSA data for JAGS run 
### 3. Format biopsy result (reclassification) data for JAGS run 




### 0 .Load libraries
library("lme4")
library("splines")


### 1. Format pt-level data for JAGS run 
#number of patients
(n <- dim(pt.data)[1])

#list of known true states
eta.true <- pt.data$true.pgg[!is.na(pt.data$true.pgg)]
#number with known true state
n_eta_known <- length(eta.true)

### removing true state observation for CV leave-one-out fold
# assign CV groups for 10 folds
set.seed(44)
my.sample <- sample(1:n_eta_known)
folds <- cut(seq(1, n_eta_known), breaks=10, labels=F)

pt.data$cvgroup<-rep(0,n)
for(group in 1:10){
  pt.data$cvgroup[my.sample[folds==group]]<-group}

#save number of patients with true state to-be-masked
(n_mask <- sum(pt.data$cvgroup==to.mask))

#save a list of subject ids, true states 
save <- as.data.frame(cbind(pt.data$subj[pt.data$cvgroup==to.mask], 
                            pt.data$true.pgg[pt.data$cvgroup==to.mask]))
names(save) <- c("subj", "true.pgg")
write.csv(save,
          paste0(location.of.generated.files,"/eta-subj-",to.mask,".csv"))

##remove true state for indicated group
pt.data$true.pgg.bin[pt.data$cvgroup==to.mask]<- NA
pt.data$true.pgg[pt.data$cvgroup==to.mask]<-NA

#reorder patients
pt.data<-pt.data[order(pt.data$true.pgg),] #order wrt binary eta still (shouldn't matter)
eta.data<-pt.data$true.pgg

#new number with known true state
(n_eta_known <- sum(!is.na(eta.data)))

#renumber patients
pt.data$subj2<-c(1:n)

#save unique ids for patients we are testing
#eta.hat.test<-pt.data$subj2[pt.data$cvgroup==to.mask]-n_eta_known

psa.data$subj2<-vector(length=dim(psa.data)[1])
for(i in 1:n){psa.data$subj2[psa.data$subj==i]<-pt.data$subj2[pt.data$subj==i]}
bx.full$subj2<-vector(length=dim(bx.full)[1])
for(i in 1:n){bx.full$subj2[bx.full$subj==i]<-pt.data$subj2[pt.data$subj==i]}




#covariate matrix for predictors of PGG
#### NEED TO ADD MAT FOR RANDOM EFFECTS MODEL
pt.data$dx.age.std <- scale(pt.data$dx.age)  ## 04/13/22 Zitong
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
subj_psa <- psa.data$subj2
sum(is.na(subj_psa))
length(unique(subj_psa))

#covariate matrix for random effects
# psa.data$years <- psa.data$days/365.25
# Z.data <- as.matrix(cbind(rep(1,n_obs_psa), psa.data$years))
# d.Z <- dim(Z.data)[2]
# data.check(condition=as.logical(sum(is.na(Z.data))==0), message="Missing time since dx in PSA data. Email Yates; she will check orginal script.")
# summary(psa.data$days)
Z.data <- as.matrix(cbind(rep(1,n_obs_psa), psa.data$time.since.dx))
d.Z <- dim(Z.data)[2]
data.check(condition=as.logical(sum(is.na(Z.data))==0), message="Missing time since dx in PSA data. Email Yates; she will check orginal script.")


#covariate matrix for fixed effects
# psa.data$dx.age.std<-scale(psa.data$dx.age.std)
# X.data <- as.matrix(cbind(psa.data$vol.std, psa.data$dx.age.std))
# d.X <- dim(X.data)[2]
# data.check(condition=as.logical(sum(is.na(X.data))==0), message="Missing volumes in PSA data. Email Yates; she will check orginal script.")
# summary(psa.data$vol.std)
# summary(psa.data$dx.age.std)
psa.data$dx.age.std<-vector(length=n_obs_psa)
for(i in 1:n){
  psa.data$dx.age.std[psa.data$subj==pt.data$subj[i]]<-pt.data$dx.age.std[i]}

#define covariate matrix
X.data <- as.matrix(cbind(psa.data$std.vol, scale(psa.data$dx.age.std)))
d.X <- dim(X.data)[2]
data.check(condition=as.logical(sum(is.na(X.data))==0), message="Missing volumes in PSA data. Email Yates; she will check orginal script.")


#lmer fit to get starting value for covariance parameter
mod_lmer<-lmer(log.psa~ std.vol + dx.age.std + (1+ time.since.dx |subj), data=psa.data)
var_vec <- apply(coef(mod_lmer)$subj, 2, var)[1:d.Z]
names(var_vec)
#make sure order of variance estimates is correct
index.intercept<-c(1:2)[names(var_vec)=="(Intercept)"]
index.time<-c(1:2)[names(var_vec)=="time.since.dx"]
var_vec <- c(var_vec[index.intercept], var_vec[index.time])


### 3. Format biopsy result (reclassification) data for JAGS run 
#number of biopsies with GS outcome
bx.full$pgg[bx.full$npc==0 & !is.na(bx.full$npc)]<-1

#subset data to intervals with biopsy occurring
pgg.data <- bx.full[bx.full$bx.here==1 & !is.na(bx.full$bx.here)
                    & !is.na(bx.full$pgg)
                    &  bx.full$time.int>0,]

n_pgg <- dim(pgg.data)[1]

#GS outcomes for biopsies
PGG <- as.numeric(pgg.data$pgg)
table(PGG)
data.check(condition=as.logical(sum(is.na(PGG))==0), message="Missing biopsy results data. Email Yates; she will check orginal script.")

#unique subj identifier
subj_pgg <- pgg.data$subj2

#covariate matrix
#covariate matrix
V.PGG.data <- as.matrix(cbind(pgg.data$bx.dt.num, pgg.data$ncs, pgg.data$mri, pgg.data$std.vol))
d.V.PGG <- dim(V.PGG.data)[2]

pgg.data2 <- pgg.data[!is.na(pgg.data$bx.time.min),]
#& !is.na(pgg.data$npc.min) & pgg.data$npc.min>=0,]
n_pgg2 <- dim(pgg.data2)[1]
PGG2 <- rep(1,n_pgg2)
subj_pgg2 <- pgg.data2$subj2
V.PGG2.data <- as.matrix(cbind(pgg.data2$bx.dt.num.min, pgg.data2$ncs.min,
                               pgg.data2$mri.min, pgg.data2$std.vol))

#combine all outcomes
n_pgg <- n_pgg + n_pgg2
PGG <- c(PGG, PGG2)
subj_pgg <- c(subj_pgg, subj_pgg2)
data.check(condition=as.logical(sum(is.na(PGG))==0), message="Missing biopsy results data. Email Yates; she will check orginal script.")

V.PGG.data <- rbind(V.PGG.data, V.PGG2.data)
#natural splines for continuous variables
bx.date.knots <- attr(ns(V.PGG.data[,1], 3), "knots")
bx.date.bknots <- attr(ns(V.PGG.data[,1], 3), "Boundary.knots")
pgg.ncs.knots <- attr(ns(V.PGG.data[,2], 2), "knots")
pgg.ncs.bknots <- attr(ns(V.PGG.data[,2], 2), "Boundary.knots")

V.PGG.data <- as.matrix(cbind( ns(V.PGG.data[,1], 3),
                               ns(V.PGG.data[,2], 2),
                               V.PGG.data[,3:d.V.PGG] ))
d.V.PGG <- dim(V.PGG.data)[2]
data.check(condition=as.logical(sum(is.na(V.PGG.data))==0), message="Missing biopsy data. Email Yates; she will check orginal script.")









