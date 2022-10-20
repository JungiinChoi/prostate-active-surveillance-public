source("~/Downloads/prostate-active-surveillance-vDaan/R-script/load_libs.R")
loc_coef <- "~/Downloads/prostate-active-surveillance-vDaan/generated-files-seq/"
load("~/Downloads/prostate-active-surveillance-vDaan/generated-files-seq/IOP-data-shaping-work-space.RData")
### Data check function
options(warn=1)
data.check <- function(condition, message){
  if(condition==FALSE){print(paste(message, "Program terminated.", sep=" "))}
  stopifnot(condition)}
source("~/Downloads/prostate-active-surveillance-vDaan/R-script/data-prep-for-jags-reform.R")

cancer_int1 <- read_csv(paste0(loc_coef, "jags-prediction-cancer_int1-2022.csv"))
cancer_int2 <- read_csv(paste0(loc_coef, "jags-prediction-cancer_int2-2022.csv"))
cancer_int3 <- read_csv(paste0(loc_coef, "jags-prediction-cancer_int3-2022.csv"))
cancer_slope1 <- read_csv(paste0(loc_coef, "jags-prediction-cancer_slope1-2022.csv"))
cancer_slope2 <- read_csv(paste0(loc_coef, "jags-prediction-cancer_slope2-2022.csv"))
cancer_slope3 <- read_csv(paste0(loc_coef, "jags-prediction-cancer_slope3-2022.csv"))
expit <- function(x){exp(x)/(1+exp(x))}
canc_int1 <- as.matrix(cancer_int1[,-1])
canc_int2 <- as.matrix(cancer_int2[,-1])
canc_int3 <- as.matrix(cancer_int3[,-1])
canc_slope1 <- as.matrix(cancer_slope1[,-1])
canc_slope2 <- as.matrix(cancer_slope2[,-1])
canc_slope3 <- as.matrix(cancer_slope3[,-1])

# canc_int1 <- apply(canc_int1,2,mean)
# canc_int2 <- apply(canc_int2,2,mean)
# canc_int3 <- apply(canc_int3,2,mean)
# canc_slope1 <- apply(canc_slope1,2,mean)
# canc_slope2 <- apply(canc_slope2,2,mean)
# canc_slope3 <- apply(canc_slope3,2,mean)

expon1 <- expit(as.vector(canc_int1) + t(modmat_cancer %*% t(canc_slope1)))
expon2 <- expit(as.vector(canc_int2) + t(modmat_cancer %*% t(canc_slope2)))
expon3 <- expit(as.vector(canc_int3) + t(modmat_cancer %*% t(canc_slope3)))

p_eta1_given_coef <- expon1
p_eta2_given_coef <- (1-expon1) * expon2
p_eta3_given_coef <- (1-expon1) * (1-expon2) *expon3
p_eta4_given_coef <- 1-(p_eta1_given_coef + p_eta2_given_coef + p_eta3_given_coef)


pgg_int1 <- read_csv(paste0(loc_coef, "jags-prediction-pgg_int1-2022.csv"))
pgg_int2 <- read_csv(paste0(loc_coef, "jags-prediction-pgg_int2-2022.csv"))
pgg_int3 <- read_csv(paste0(loc_coef, "jags-prediction-pgg_int3-2022.csv"))
pgg_slope1 <- read_csv(paste0(loc_coef, "jags-prediction-pgg_slope1-2022.csv"))
pgg_slope2 <- read_csv(paste0(loc_coef, "jags-prediction-pgg_slope2-2022.csv"))
pgg_slope3 <- read_csv(paste0(loc_coef, "jags-prediction-pgg_slope3-2022.csv"))

pg_int1 <- as.matrix(pgg_int1[,-1])
pg_int2 <- as.matrix(pgg_int2[,-1])
pg_int3 <- as.matrix(pgg_int3[,-1])
pg_slope1 <- as.matrix(pgg_slope1[,-1])
pg_slope2 <- as.matrix(pgg_slope2[,-1])
pg_slope3 <- as.matrix(pgg_slope3[,-1])

## update covariate matrix for biopsy model with eta predictions
K = 4
eta.data<-pt.data$true.pgg
etahat<-read.csv(paste0(loc_coef, "/jags-prediction-cancer_state-2022.csv"))
etahat<-as.matrix(etahat[,2:dim(etahat)[2]])
(n <- dim(pt.data)[1])
(B<-dim(etahat)[1]) #length of saved posterior chain
#mean eta predictions for patients with unknown cancer state
pred_eta <- matrix(nrow=(n-npat_cancer_known), ncol=K)
for(i in 1:(n-npat_cancer_known)){for(k in 1:K){
  pred_eta[i,k] <- sum(etahat[,i]==k)/B }}

V.PGG.eta<-matrix(nrow=npat_pgg, ncol=(K-1))
for(j in 1:npat_pgg){
  if(pgg_patient_index_map[j] <= npat_cancer_known){
    for(k in 1:(K-1)){
      V.PGG.eta[j,k] <- as.numeric(eta.data[pgg_patient_index_map[j]]>k)}}
  else{
    for(k in 1:(K-1)){
      V.PGG.eta[j,k] <- 1-sum(pred_eta[(pgg_patient_index_map[j]-npat_cancer_known),1:k]) } } }

modmat_pg <- cbind(modmat_pgg, V.PGG.eta)
expon_pg1 <- expit(as.vector(pg_int1) + t(modmat_pg %*% t(pg_slope1)))
expon_pg2 <- expit(as.vector(pg_int2) + t(modmat_pg %*% t(pg_slope2)))
expon_pg3 <- expit(as.vector(pg_int3) + t(modmat_pg %*% t(pg_slope3)))

p_pg1_given_coef <- expon_pg1
p_pg2_given_coef <- (1-expon_pg1) * expon_pg2
p_pg3_given_coef <- (1-expon_pg1) * (1-expon_pg2) *expon_pg3
p_pg4_given_coef <- 1-(p_pg1_given_coef + p_pg2_given_coef + p_pg3_given_coef)

fixef_psa <- read_csv(paste0(loc_coef, "/jags-prediction-fixef_coefficient-2022.csv"))
ranef_int <- read_csv(paste0(loc_coef, "/jags-prediction-ranef_intercept-2022.csv"))
ranef_slope <- read_csv(paste0(loc_coef, "/jags-prediction-ranef_slope-2022.csv"))

ff_psa <- as.matrix(fixef_psa[,-1])
rf_psa_int <- as.matrix(ranef_int[,-1])
rf_psa_slope <- as.matrix(ranef_slope[,-1])

lmm_mean <- modmat_fixef_psa %*% t(ff_psa) 
var_ranef_int <- read_csv(paste0(loc_coef, "/jags-prediction-ranef_var_intercept-2022.csv"))
var_ranef_slope <- read_csv(paste0(loc_coef, "/jags-prediction-ranef_var_slope-2022.csv"))
cov_ranef <- read_csv(paste0(loc_coef, "/jags-prediction-ranef_cov-2022.csv"))
var_resid <- read_csv(paste0(loc_coef, "/jags-prediction-resid_var_psa-2022.csv"))
var_rf_int <- as.matrix(var_ranef_int[,-1])
var_rf_slope <- as.matrix(var_ranef_slope[,-1])
cov_rf <- as.matrix(cov_ranef[,-1])
var_e <- as.matrix(var_resid[,-1])
(B <- nrow(cov_rf))
# vcov_rf <- list()
# psa_marginal_var <- list()
# for(i in 1:B){
#   vcov_rf[[i]] <- matrix(c(var_rf_int[i], cov_rf[i], cov_rf[i],var_rf_slope[i]), nrow = 2, byrow = T)
#   psa_marginal_var[[i]] <- modmat_ranef_psa %*% vcov_rf[[i]] %*% t(modmat_ranef_psa) + diag(var_e[i],nrow(modmat_ranef_psa))
# }
sd_rf_int <- apply(var_rf_int, 2, mean)
sd_rf_slope <- apply(var_rf_slope, 2, mean)
cov_rf <- apply(cov_rf, 2, mean)
sd_e <- apply(var_e, 2, mean)
vcov_rf <- matrix(c(sd_rf_int^2, cov_rf, cov_rf,sd_rf_slope^2), nrow = 2, byrow = T)
#psa_marginal_var <- modmat_ranef_psa %*% vcov_rf %*% t(modmat_ranef_psa) + diag(sd_e^2,nrow(modmat_ranef_psa))

rf_eta12 <- cbind(rf_psa_int[,1], rf_psa_slope[,1])
rf_eta34 <- cbind(rf_psa_int[,2], rf_psa_slope[,2])
psa_mean_eta12 <- modmat_fixef_psa %*% t(ff_psa)  + modmat_ranef_psa %*%t(rf_eta12)
psa_mean_eta34 <- modmat_fixef_psa %*% t(ff_psa)  + modmat_ranef_psa %*%t(rf_eta34)

psa_obs <- psa.data$log.psa

for(i in 1:)
lik_psa_eta12 <- (2*pi)^()

