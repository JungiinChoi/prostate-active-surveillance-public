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
cancer_intercept  <- matrix(c(0.4, 2,4))
cancer_coef <- matrix(c(0.5, 0.8))
## mu_k for mean of b
mu_eta0 <- matrix(c(1.3, 0.02))
mu_eta1 <- matrix(c(1.5,0.1))
## beta, Sgima, alpha, gamma and sigma^2
Beta <- matrix(c(0.5, 0.01)) #fixef coefficients for lme
Sigma <- matrix(c(0.5, 0.01, 0.01, 0.01), nrow = 2) #ranef vcov for lme
sigma2 <- 0.3 #variance of Y for lme
alpha <- matrix(c(3.25, 4.58, 5.84))
gamma <- matrix(c(0.6, -2, 0.9, 2, -8, 0.5, -0.3, 1.9, 1.3, 0.5)) # coefs for prop odds model

# simulate eta----
simulate_cancer_state <- function(intercept, coef, predictor) {
  pred_effect <- predictor %*% coef
  cum_prob <- lapply(
    1:3, function(k) {
      logit_cum_prob <- intercept[k] - pred_effect
      return(1 / (1 + exp(-logit_cum_prob)))
    }
  )
  prob <- cbind(
    cum_prob[[1]],
    cum_prob[[2]] - cum_prob[[1]],
    cum_prob[[3]] - cum_prob[[2]],
    1 - cum_prob[[3]]
  )
  cancer_state <- sapply(
    1:dim(predictor)[1], 
    function(i) sample(1:4, size=1, prob=prob[i, ])
  )
  return(cancer_state)
}

set.seed(2021)
eta.sim <- simulate_cancer_state(cancer_intercept, cancer_coef, V.ETA.data)

## create missings for eta
set.seed(2021)
obs.eta <- sample(1:n, length(eta.data))
eta.data.sim <- eta.sim[obs.eta]

# simulate e, b and Y ------ 
generate_rand_coef <- function(cancer_state, mean_0, mean_1, cov) {
  rand_coef_mean <- if (cancer_state <= 2) mean_0 else mean_1
  rand_coef <- rmvnorm(1, rand_coef_mean, cov)
  return(rand_coef)
}

simulate_psa <- function(
    cancer_state, fixed_eff_pred, rand_eff_pred, fixed_coef, rand_coef, 
    resid_var, psa_patient_index_map
  ) {
  # Predictors are in the long-format, with both patient and time along the column.
  # Incidentally, the fixed effect predictors are *not* time-dependent. And the
  # intercept is treated as a random effect, as opposed to what's documented in
  # the manuscript.
  n_psa_meas <- length(psa_patient_index_map)
  psa_mean <- sapply(
    1:n_psa_meas, function(j) {
      rand_coef_j <- rand_coef[[psa_patient_index_map[j]]] 
      fixed_coef %*% fixed_eff_pred[j, ] + rand_coef_j %*% rand_eff_pred[j, ]
    }
  )
  psa <- psa_mean + rnorm(n_psa_meas, sd = sqrt(resid_var))
  return(psa)
}

set.seed(2021)
rand_coef <- lapply(
  1:n, function(i) generate_rand_coef(eta.sim[i], mu_eta0, mu_eta1, Sigma)
)
resid_var <- 0.5 
psa <- simulate_psa(eta.sim, X.data, Z.data, t(Beta), rand_coef, resid_var, subj_psa)

# simulate R(PGG) ---
simulate_pgg_state <- function(
  cancer_state, biopsy_pred, intercept, cov_eff_coef, cancer_eff_coef, 
  pgg_patient_index_map
) {
  pred_effect <- biopsy_pred %*% cov_eff_coef
  for (k in 1:3) {
    grade_above_k <- as.numeric(cancer_state[pgg_patient_index_map] > k)
    pred_effect <- pred_effect + grade_above_k * cancer_eff_coef[k]
  }
  cum_prob <- lapply(
    1:3, function(k) {
      logit_cum_prob <- intercept[k] - pred_effect
      return(1 / (1 + exp(-logit_cum_prob)))
    }
  )
  pgg_prob <- cbind(
    cum_prob[[1]],
    cum_prob[[2]] - cum_prob[[1]],
    cum_prob[[3]] - cum_prob[[2]],
    1 - cum_prob[[3]]
  )
  pgg_state <- sapply(
    1:dim(biopsy_pred)[1],
    function(i) sample(1:4, size=1, prob=pgg_prob[i, ])
  )
  return(pgg_prob)
}

pgg_prob <- simulate_pgg_state(
  eta.sim, V.PGG.data, alpha, gamma[1:7], gamma[8:10], subj_pgg
)

PGG.sim <- list()
pgg_prob_prev <- list()
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
  pgg_prob_i <- NULL
  for(j in 1:length(inx)){
    tmp <- sample(1:4, 1, replace=TRUE, prob= c(p.pgg1[j], p.pgg2[j], p.pgg3[j], p.pgg4[j]) )
    tmp.PGG <- rbind(tmp.PGG, tmp)
    pgg_prob_i <- rbind(pgg_prob_i, c(p.pgg1[j], p.pgg2[j], p.pgg3[j], p.pgg4[j]))
  }
  PGG.sim[[i]] <- matrix(tmp.PGG, ncol = 1)
  pgg_prob_prev[[i]] <- pgg_prob_i
}
PGG.new <- do.call("rbind", PGG.sim)
pgg_prob_prev <- do.call(rbind, pgg_prob_prev)

stopifnot(all(abs(pgg_prob - pgg_prob_prev) < 1e-8)) # Check that the two outputs coincide

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
plot(jags_rhoint$V1, type = "l")
