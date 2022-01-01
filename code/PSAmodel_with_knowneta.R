# Proportional odds and mixed modelfor patients 
#with known eta (only people with prostectmy)
#Set prior to those estimates and set not too strong variance, compare jags code with these 
#and compare jags using simulation with missing eta ~ 5% 
#See if estimation agrees with each other

rm(list=ls())
source('R-scripts/load_libs.R')
base.location <- "/Users/zitongwang/Downloads/prostate-active-surveillance-vDaan/" #"/users/rcoley/jhas-epic/psa/psa-new/" #"/Users/ryc/Documents/inhealth/prediction-model/automated/for-TIC/"
setwd(base.location)
location.of.data <- paste0(base.location, "data")
location.of.r.scripts <- paste0(base.location, "R-scripts")
location.of.generated.files <- paste0(base.location, "generated-files-etaknown-model")

name.of.pt.data <- "Demographic_6.15.csv" #"demographics with physician info.2015.csv"
name.of.bx.data <- "Biopsy_6.15.csv" #"Biopsy data_2015.csv"
name.of.psa.data <- "PSA_6.15.csv" #"PSA.2015.csv"
name.of.tx.data <- "Treatment_6.15.csv" #"treatment_2015.csv"

date.pull<-as.numeric(Sys.Date())
source(paste(location.of.r.scripts,"data-load-check-and-shaping.R",sep="/"))
source(paste(location.of.r.scripts,"data-prep-for-jags.R",sep="/"))
source(paste(location.of.r.scripts,"argument-prep-for-jags.R",sep="/"))
source(paste(location.of.r.scripts,"JAGS-prediction-model-etaknown.R",sep="/"))

eta.info <- pt.data %>% dplyr::select(clinical_PTnum,subj, true.pgg)
psa.data.eta <- psa.data %>% left_join(eta.info)
psa.data.eta <- psa.data.eta %>% filter(!is.na(true.pgg))
mod_lmer<-lmer(log.psa~ std.vol + dx.age.std + (1+ time.since.dx |subj), data=psa.data.eta)

# use only known eta patients
pt.data.eta <- pt.data %>% filter(!is.na(true.pgg))
n <- dim(pt.data.eta)[1]
eta.data <- pt.data.eta$true.pgg[!is.na(pt.data.eta$true.pgg)]
n_eta_known <- length(eta.data)
pt.data.eta$dx.age.std <- scale(pt.data.eta$dx.age)
V.ETA.data <- as.matrix(cbind(pt.data.eta$dx.age.std, pt.data.eta$lr.vol))
d.V.ETA<-dim(V.ETA.data)[2]

#psa
n_obs_psa <- dim(psa.data.eta)[1]
Y <- psa.data.eta$log.psa
data.check(condition=as.logical(sum(is.na(Y))==0), message="Missing PSA values. Email Yates; she will check orginal script.")
#list of patients
subj_psa <- psa.data.eta$subj

#covariate matrix for random effects
Z.data <- as.matrix(cbind(rep(1,n_obs_psa), psa.data.eta$time.since.dx))
d.Z <- dim(Z.data)[2]
data.check(condition=as.logical(sum(is.na(Z.data))==0), message="Missing time since dx in PSA data. Email Yates; she will check orginal script.")
psa.data.eta$dx.age.std<-vector(length=n_obs_psa)
for(i in 1:n){
  psa.data.eta$dx.age.std[psa.data.eta$subj==
                            pt.data.eta$subj[i]]<-pt.data.eta$dx.age.std[i]}
X.data <- as.matrix(cbind(psa.data.eta$std.vol, scale(psa.data.eta$dx.age.std)))
d.X <- dim(X.data)[2]
data.check(condition=as.logical(sum(is.na(X.data))==0), message="Missing volumes in PSA data. Email Yates; she will check orginal script.")
var_vec <- apply(coef(mod_lmer)$subj, 2, var)[1:d.Z]
names(var_vec)
index.intercept<-c(1:2)[names(var_vec)=="(Intercept)"]
index.time<-c(1:2)[names(var_vec)=="time.since.dx"]
var_vec <- c(var_vec[index.intercept], var_vec[index.time])
#biospy result
bx.full.eta <- bx.full %>% left_join(eta.info)
bx.full.eta <- bx.full.eta %>% filter(!is.na(true.pgg))
bx.full.eta$pgg[bx.full.eta$npc==0 & !is.na(bx.full.eta$npc)]<-1
pgg.data <- bx.full.eta[bx.full.eta$bx.here==1 & !is.na(bx.full.eta$bx.here)
                    & !is.na(bx.full.eta$pgg)
                    &  bx.full.eta$time.int>0,]
n_pgg <- dim(pgg.data)[1]
PGG <- as.numeric(pgg.data$pgg)
subj_pgg <- pgg.data$subj
V.PGG.data <- as.matrix(cbind(pgg.data$bx.dt.num, pgg.data$ncs, pgg.data$mri, pgg.data$std.vol))
d.V.PGG <- dim(V.PGG.data)[2]
pgg.data2 <- pgg.data[!is.na(pgg.data$bx.time.min),]
n_pgg2 <- dim(pgg.data2)[1]
PGG2 <- rep(1,n_pgg2)
subj_pgg2 <- pgg.data2$subj
V.PGG2.data <- as.matrix(cbind(pgg.data2$bx.dt.num.min, pgg.data2$ncs.min,
                               pgg.data2$mri.min, pgg.data2$std.vol))

#combine all outcomes
n_pgg <- n_pgg + n_pgg2
PGG <- c(PGG, PGG2)
subj_pgg <- c(subj_pgg, subj_pgg2)
data.check(condition=as.logical(sum(is.na(PGG))==0), message="Missing biopsy results data. Email Yates; she will check orginal script.")

V.PGG.data <- rbind(V.PGG.data, V.PGG2.data)
bx.date.knots <- attr(ns(V.PGG.data[,1], 3), "knots")
bx.date.bknots <- attr(ns(V.PGG.data[,1], 3), "Boundary.knots")
pgg.ncs.knots <- attr(ns(V.PGG.data[,2], 2), "knots")
pgg.ncs.bknots <- attr(ns(V.PGG.data[,2], 2), "Boundary.knots")

V.PGG.data <- as.matrix(cbind( ns(V.PGG.data[,1], 3),
                               ns(V.PGG.data[,2], 2),
                               V.PGG.data[,3:d.V.PGG] ))
d.V.PGG <- dim(V.PGG.data)[2]
data.check(condition=as.logical(sum(is.na(V.PGG.data))==0), message="Missing biopsy data. Email Yates; she will check orginal script.")

seed = 2021
set.seed(seed)
jags_data<-list(K=K, K.bin=K.bin, n=n,
                eta.data=eta.data, n_eta_known=n_eta_known,
                V.ETA=V.ETA.data, d.V.ETA=d.V.ETA,
                
                n_obs_psa=n_obs_psa, Y=Y, subj_psa=subj_psa,
                Z=Z.data, X=X.data, d.Z=d.Z, d.X=d.X, I_d.Z=diag(d.Z),
                
                PGG=PGG, n_pgg=n_pgg, subj_pgg=subj_pgg,
                V.PGG=V.PGG.data, d.V.PGG=d.V.PGG#,
                
) #,
outj<-jags(jags_data, inits=inits, parameters.to.save=params,
           model.file=paste(location.of.r.scripts, "JAGS-prediction-model-etaknown.txt", sep="/"),
           n.thin=n.thin, n.chains=n.chains, n.burnin=n.burnin, n.iter=n.iter)

out<-outj$BUGSoutput
for(j in 1:length(out$sims.list)){
  write.csv(out$sims.list[[j]],
            paste(location.of.generated.files, "/jags-prediction-", names(out$sims.list)[j],"-",seed,".csv",sep=""))
}

library(MASS)
pom.PGG.data <- as.data.frame(V.PGG.data)
pgg.data.all <- rbind(pgg.data, pgg.data2)
pom.PGG.data$eta_m2 <- ifelse(pgg.data.all$true.pgg-2 >=0 , 1, 0)
pom.PGG.data$eta_m3 <- ifelse(pgg.data.all$true.pgg-3 >=0 , 1, 0)
pom.PGG.data$eta_m4 <- ifelse(pgg.data.all$true.pgg-4 >=0 , 1, 0)
pom.PGG.data <- as.matrix(pom.PGG.data)
mod_pom <- polr(factor(PGG) ~ pom.PGG.data)
summary(mod_pom)
summary(mod_lmer)

get_stats <- function(x){a <- quantile(x, c(0.025, 0.975)); b<-mean(x); return(c(a[1], b, a[2]))}
jags_beta <- read_csv(paste0(location.of.generated.files,"/jags-prediction-beta-2021.csv"))
apply(jags_beta, 2, get_stats)
jags_alpha <- read_csv(paste0(location.of.generated.files,"/jags-prediction-alpha-2021.csv"))
apply(jags_alpha, 2, get_stats)
jags_gamma<- read_csv(paste0(location.of.generated.files,"/jags-prediction-gamma.PGG-2021.csv"))
apply(jags_gamma, 2, get_stats)


## jags with simulated data -----
jags_rho_int<- read_csv(paste0(location.of.generated.files,"/jags-prediction-rho_int-2021.csv"))
jags_rho_coef <- read_csv(paste0(location.of.generated.files,"/jags-prediction-rho_coef-2021.csv"))
jags_mu_int <- read_csv(paste0(location.of.generated.files,"/jags-prediction-mu_int-2021.csv"))
jags_mu_slp <- read_csv(paste0(location.of.generated.files,"/jags-prediction-mu_slope-2021.csv"))
jags_beta <- read_csv(paste0(location.of.generated.files,"/jags-prediction-beta-2021.csv"))
jags_D_int <- read_csv(paste0(location.of.generated.files,"/jags-prediction-sigma_int-2021.csv"))
jags_D_slp <- read_csv(paste0(location.of.generated.files,"/jags-prediction-sigma_slope-2021.csv"))
jags_D_cov <- read_csv(paste0(location.of.generated.files,"/jags-prediction-cov_int_slope-2021.csv"))
jags_e <- read_csv(paste0(location.of.generated.files,"/jags-prediction-sigma_res-2021.csv"))
jags_alpha <- read_csv(paste0(location.of.generated.files,"/jags-prediction-alpha-2021.csv"))
jags_gamma <- read_csv(paste0(location.of.generated.files,"/jags-prediction-gamma.PGG-2021.csv"))


apply(jags_rho_int, 2, get_stats)
apply(jags_rho_coef, 2, get_stats)
apply(jags_mu_int, 2, get_stats)
apply(jags_mu_slp, 2, get_stats)
apply(jags_beta, 2, get_stats)
apply(jags_D_int, 2, get_stats)
apply(jags_D_slp, 2, get_stats)
apply(jags_D_cov, 2, get_stats)
apply(jags_e, 2, get_stats)
apply(jags_alpha, 2, get_stats)
apply(jags_gamma, 2, get_stats)


location.of.generated.files <- paste0(base.location, "generated-files")

## prespecify rho intercept and rho slope
rho_int  <- matrix(c( -0.16, 1.62, 3.28))
rho_coef <- matrix(c(0.48, 0.06))
## mu_k for mean of b
mu_eta0 <- matrix(c(1.47, 0.05))
mu_eta1 <- matrix(c(1.54,0.11))
## beta, Sgima, alpha, gamma and sigma^2
Beta <- matrix(c(0.2, 0.09)) #fixef coefficients for lme
Sigma <- matrix(c(0.45, 0.01, 0.01, 0.1), nrow = 2) #ranef vcov for lme
sigma2 <- 0.27 #variance of Y for lme
alpha <- matrix(c(1.77, 3.17, 4.55))
gamma <- matrix(c(0.85, -1.39, 1.94, 0.29, 1.28, 0.33, -0.23, 0.89, 1.19, 0.28)) # coefs for prop odds model

# simulate eta------
eta.linpred <- V.ETA.data  %*% rho_coef
logit.cum.e1 <- rho_int[1] - eta.linpred  
logit.cum.e2 <- rho_int[2] - eta.linpred  
logit.cum.e3 <- rho_int[3] - eta.linpred  

cum.e1 <- exp(logit.cum.e1)/(1+exp(logit.cum.e1))
cum.e2 <- exp(logit.cum.e2)/(1+exp(logit.cum.e2))
cum.e3 <- exp(logit.cum.e3)/(1+exp(logit.cum.e3))

p.eta4 <- 1-cum.e3
p.eta3 <- cum.e3 - cum.e2
p.eta2 <- cum.e2 - cum.e1
p.eta1 <- cum.e1
set.seed(seed)
tmp.eta <- NULL
for(l in 1:n){ 
  tmp <- sample(1:4, 1, replace=TRUE, prob= c(p.eta1[l], p.eta2[l], p.eta3[l], p.eta4[l]) )
  tmp.eta <- rbind(tmp.eta, tmp)
}
eta.sim <- matrix(tmp.eta, ncol = 1)

## (removed) create missings for eta
eta.data.sim <- c(eta.sim)

# simulate e, b and Y ------
set.seed(2021)
resid_var <- 0.5   #110121
b <- NULL
for(i in 1:n){
  if(eta.sim[i] <= 2){ ## use eta.bin here rather than eta
    b.tmp = rmvnorm(1, mu_eta0, Sigma)
  }else{
    b.tmp = rmvnorm(1, mu_eta1, Sigma)
  }
  b <- rbind(b, b.tmp)
}

Y.sim <- list()
for(i in 1:length(unique(subj_psa))){
  j <- unique(subj_psa)[i]
  inx <- which(subj_psa == j)
  etai <- eta.sim[j]
  Xi <- X.data[inx, ]
  bi <- matrix(b[j,])
  ei <- resid_var
  Zi <- Z.data[inx, ]
  Yi.mean = Xi %*% Beta + Zi %*% bi
  Yi = rmvnorm(1, Yi.mean, diag(ei, nrow = length(inx)))
  Y.sim[[i]] <- matrix(Yi, ncol = 1)
}
Y.new <- do.call("rbind", Y.sim)

# simulate R(PGG) ---
PGG.sim <- list()
for(i in 1:length(unique(subj_pgg))){
  j <- unique(subj_pgg)[i]
  inx <- which(subj_pgg == j)
  etai <- eta.sim[j]
  Vi.pgg <- matrix(V.PGG.data[inx,], ncol = d.V.PGG)
  step.etasub2 <- ifelse(etai - 2 >=0, 1, 0)
  step.etasub3 <- ifelse(etai - 3 >=0, 1, 0)
  step.etasub4 <- ifelse(etai - 4 >=0, 1, 0)
  tmp.Vi.pgg <- cbind(Vi.pgg, 
                      matrix(rep(step.etasub2, nrow(Vi.pgg))),
                      matrix(rep(step.etasub3, nrow(Vi.pgg))),
                      matrix(rep(step.etasub4, nrow(Vi.pgg))))
  lin.pred <- tmp.Vi.pgg %*% gamma
  logit.cum.t1 <- alpha[1] - lin.pred  #logit(R <= 1) = logit(R = 1)
  logit.cum.t2 <- alpha[2] - lin.pred  #logit(R <= 2) = logit(R=1 | 2)
  logit.cum.t3 <- alpha[3] - lin.pred  #logit(R <= 3) = logit(R=1 | 2 | 3) 
  
  cum.t1 <- exp(logit.cum.t1)/(1+exp(logit.cum.t1))
  cum.t2 <- exp(logit.cum.t2)/(1+exp(logit.cum.t2))
  cum.t3 <- exp(logit.cum.t3)/(1+exp(logit.cum.t3))
  
  p.pgg4 <- 1-cum.t3
  p.pgg3 <- cum.t3 - cum.t2
  p.pgg2 <- cum.t2 - cum.t1
  p.pgg1 <- cum.t1
  
  tmp.PGG <- NULL
  for(j in 1:length(inx)){
    tmp <- sample(1:4, 1, replace=TRUE, prob= c(p.pgg1[j], p.pgg2[j], p.pgg3[j], p.pgg4[j]) )
    tmp.PGG <- rbind(tmp.PGG, tmp)
  }
  PGG.sim[[i]] <- matrix(tmp.PGG, ncol = 1)
}
PGG.new <- do.call("rbind", PGG.sim)


### run model -----
jags_data_sim<-list(K=K, K.bin=K.bin, n=n,
                eta.data=eta.data.sim, n_eta_known=n_eta_known,
                V.ETA=V.ETA.data, d.V.ETA=d.V.ETA,
                
                n_obs_psa=n_obs_psa, Y=Y.new[,1], subj_psa=subj_psa,
                Z=Z.data, X=X.data, d.Z=d.Z, d.X=d.X, I_d.Z=diag(d.Z),
                
                PGG=PGG.new[,1], n_pgg=n_pgg, subj_pgg=subj_pgg,
                V.PGG=V.PGG.data, d.V.PGG=d.V.PGG#,
)
#do.one(2021)
seed = 2021
set.seed(seed)
outj_sim<-jags(jags_data_sim, inits=inits, parameters.to.save=params,
           model.file=paste(location.of.r.scripts, "JAGS-prediction-model-etaknown.txt", sep="/"),
           n.thin=n.thin, n.chains=n.chains, n.burnin=n.burnin, n.iter=n.iter)
outj_sim<-outj_sim$BUGSoutput
for(j in 1:length(outj_sim$sims.list)){
  write.csv(outj_sim$sims.list[[j]],
            paste(location.of.generated.files, "/jags-prediction-", names(outj_sim$sims.list)[j],"-",seed,".csv",sep=""))
}
 
jags_beta <- read_csv(paste0(location.of.generated.files,"/jags-prediction-beta-2021.csv"))
apply(jags_beta, 2, get_stats)
t(Beta)
jags_alpha <- read_csv(paste0(location.of.generated.files,"/jags-prediction-alpha-2021.csv"))
apply(jags_alpha, 2, get_stats)
t(alpha)
jags_gamma<- read_csv(paste0(location.of.generated.files,"/jags-prediction-gamma.PGG-2021.csv"))
apply(jags_gamma, 2, get_stats)
t(gamma)
