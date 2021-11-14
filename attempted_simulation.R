rm(list=ls())
source('R-scripts/load_libs.R')
location.of.data <- "data"
location.of.r.scripts <- "R-scripts"
location.of.generated.files <- "generated-files-simulation-strongprior"

name.of.pt.data <- "Demographic_6.15.csv" #"demographics with physician info.2015.csv"
name.of.bx.data <- "Biopsy_6.15.csv" #"Biopsy data_2015.csv"
name.of.psa.data <- "PSA_6.15.csv" #"PSA.2015.csv"
name.of.tx.data <- "Treatment_6.15.csv" #"treatment_2015.csv"

# pt.data.original <- read.csv(paste0(location.of.data, "/",name.of.pt.data))
# bx.data.original <- read.csv(paste0(location.of.data, "/",name.of.bx.data))
# psa.data.original <- read.csv(paste0(location.of.data, "/",name.of.psa.data))
# tx.data.original <- read.csv(paste0(location.of.data, "/",name.of.tx.data))
# 
# ## bootstrap by patient (use psa.data patients) # which(!(psa.data.original$clinical_PTnum %in% pt.data.original$clinical_PTnum))
# ptnum <- unique(psa.data.original$clinical_PTnum)
# set.seed(2021)
# boot_ptnum <- sample(ptnum, 2*length(ptnum), replace = T)
# new_PTnum <- 1:length(boot_ptnum)
# pt.data <- NULL
# for(i in 1:length(boot_ptnum)){
#   tmp <- pt.data.original[pt.data.original$clinical_PTnum == boot_ptnum[i],]
#   if(nrow(tmp) != 0){tmp$new_PTnum <- new_PTnum[i]}
#   pt.data <- rbind(pt.data, tmp)
# }
# 
# 
# ## match other datasets
# bx.data <- NULL
# for(i in 1:length(boot_ptnum)){
#   tmp <- bx.data.original[which(bx.data.original$clinical_PTnum == boot_ptnum[i]),]
#   if(nrow(tmp) != 0){tmp$new_PTnum <- new_PTnum[i]}
#   bx.data <- rbind(bx.data, tmp)
# }
# 
# psa.data <- NULL
# for(i in 1:length(boot_ptnum)){
#   tmp <- psa.data.original[which(psa.data.original$clinical_PTnum == boot_ptnum[i]),]
#   if(nrow(tmp) != 0){tmp$new_PTnum <- new_PTnum[i]}
#   psa.data <- rbind(psa.data, tmp)
# }
# 
# tx.data <- NULL
# for(i in 1:length(boot_ptnum)){
#   tmp <- tx.data.original[which(tx.data.original$clinical_PTnum == boot_ptnum[i]),]
#   if(nrow(tmp) != 0){tmp$new_PTnum <- new_PTnum[i]}
#   tx.data <- rbind(tx.data, tmp)
# }
# 
# psa.data <- psa.data %>% dplyr::rename(clinical_PTnum_o = clinical_PTnum,
#                                        clinical_PTnum = new_PTnum)
# pt.data <- pt.data %>% dplyr::rename(clinical_PTnum_o = clinical_PTnum,
#                                        clinical_PTnum = new_PTnum)
# bx.data <- bx.data %>% dplyr::rename(clinical_PTnum_o = clinical_PTnum,
#                                        clinical_PTnum = new_PTnum)
# tx.data <- tx.data %>% dplyr::rename(clinical_PTnum_o = clinical_PTnum,
#                                        clinical_PTnum = new_PTnum)
# write.csv(psa.data, "data/PSA_boot.csv", row.names = F)
# write.csv(pt.data, "data/Demographic_boot.csv", row.names = F)
# write.csv(bx.data, "data/Biopsy_boot.csv", row.names = F)
# write.csv(tx.data, "data/Treatment_boot.csv", row.names = F)

date.pull<-as.numeric(Sys.Date())
source(paste(location.of.r.scripts,"data-load-check-and-shaping.R",sep="/"))
source(paste(location.of.r.scripts,"data-prep-for-jags.R",sep="/"))
source(paste(location.of.r.scripts,"argument-prep-for-jags.R",sep="/"))
source(paste(location.of.r.scripts,"JAGS-prediction-model.R",sep="/"))

# set param values ----
## prespecify rho intercept and rho slope
rho_int  <- matrix(c(0.4, 2,4))
rho_coef <- matrix(c(0.5, 0.8))
## mu_k for mean of b
mu_eta0 <- matrix(c(1.3, 0.02))
mu_eta1 <- matrix(c(1.5,0.1))
## beta, Sgima, alpha, gamma and sigma^2
Beta <- matrix(c(0.5, 0.01)) #fixef coefficients for lme
Sigma <- matrix(c(0.5, 0.01, 0.01, 0.01), nrow = 2) #ranef vcov for lme
sigma2 <- 0.3 #variance of Y for lme
alpha <- matrix(c(3.25, 4.58, 5.84))
gamma <- matrix(c(0.6, -2, 0.9, 2, -8, 0.5, -0.3, 1.9, 1.3, 0.5)) # coefs for prop odds model

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
set.seed(2021)
tmp.eta <- NULL
for(l in 1:n){ 
  tmp <- sample(1:4, 1, replace=TRUE, prob= c(p.eta1[l], p.eta2[l], p.eta3[l], p.eta4[l]) )
  tmp.eta <- rbind(tmp.eta, tmp)
}
eta.sim <- matrix(tmp.eta, ncol = 1)

## create missings for eta
set.seed(2021)
obs.eta <- sample(1:n, length(eta.data))
eta.data.sim <- eta.sim[obs.eta]

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
jags_data<-list(K=K, K.bin=K.bin, n=n,
                eta.data=eta.data.sim, n_eta_known=n_eta_known,
                V.ETA=V.ETA.data, d.V.ETA=d.V.ETA,
                
                n_obs_psa=n_obs_psa, Y=Y.new[,1], subj_psa=subj_psa,
                Z=Z.data, X=X.data, d.Z=d.Z, d.X=d.X, I_d.Z=diag(d.Z),
                
                PGG=PGG.new[,1], n_pgg=n_pgg, subj_pgg=subj_pgg,
                V.PGG=V.PGG.data, d.V.PGG=d.V.PGG
                
)
#do.one(2021)
seed = 2021
set.seed(seed)
outj<-jags(jags_data, inits=inits, parameters.to.save=params,
           model.file=paste(location.of.r.scripts, "JAGS-prediction-model-simp.txt", sep="/"),
           n.thin=n.thin, n.chains=n.chains, n.burnin=n.burnin, n.iter=n.iter)
out<-outj$BUGSoutput
for(j in 1:length(out$sims.list)){
  write.csv(out$sims.list[[j]],
            paste(location.of.generated.files, "/jags-prediction-", names(out$sims.list)[j],"-",seed,".csv",sep=""))
}
#save(out, file = "generated-files-with-bootstrapped-data/jags_output_sim.RData")

## compare results
get_stats <- function(x){a <- quantile(x, c(0.025, 0.975)); b<-mean(x); return(c(a[1], b, a[2]))}
jags_beta <- read_csv("generated-files/jags-prediction-beta-2021.csv")
apply(jags_beta, 2, get_stats)
t(Beta)

jags_gamma <- read_csv(paste0(location.of.generated.files,"/jags-prediction-gamma.PGG-2021.csv"))
apply(jags_gamma, 2, get_stats)
t(gamma)

jags_alpha <-read_csv(paste0(location.of.generated.files,"/jags-prediction-alpha-2021.csv"))
apply(jags_alpha, 2, get_stats)
t(alpha)

jags_rhoint <- read_csv(paste0(location.of.generated.files,"/jags-prediction-rho_int-2021.csv")
jags_rhoint_1000 <- read_csv("generated-files-old/jags-prediction-rho_int-1000.csv")
apply(jags_rhoint, 2, get_stats)
apply(jags_rhoint_1000, 2, get_stats)
t(rho_int)
plot(jags_rhoint$V1, type = "l")
