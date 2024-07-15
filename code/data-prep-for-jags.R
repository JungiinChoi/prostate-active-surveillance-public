########### Data processing for a hierarchical model using J sets of synthetic data. 



# Yates Coley
# rycoley@gmail.com
# 2019 January 08, modified 2019 January 30
# This script will take tidied/shaped AS data and perform additional manipulation for the JAGS model
# I've put a lot of data checks in here as fail-safes, but there shouldn't be any missing values after running the previous script.
# This updated from data-prep-for-jags.R to accomodated K-fold CV 

#### WORKFLOW
### 1. Format patients' demopgrachic  data for JAGS run
#      Include removing true state observation for CV leave-one-out fold
### 2. Format PSA data for JAGS run 
### 3. Format biopsy result (reclassification) data for JAGS run 


### 0 .Load libraries --------
library("lme4")
library("splines")
library("readr")
library("dplyr")
library("tidyr")
library("readxl")

# load data files from different institutions

# number of institutions 

dx_list <- bx_list <- mri_list <- psa_list <- list(length = J)

for (i in 1:J){
  dx_list[[i]] <- read.csv(paste0(location.of.data,"/dx_",i,".csv"))
  bx_list[[i]] <- read.csv(paste0(location.of.data,"/bx_",i,".csv"))
  psa_list[[i]] <- read.csv(paste0(location.of.data,"/psa_",i,".csv"))
  mri_list[[i]] <- read.csv(paste0(location.of.data,"/mri_",i,".csv"))
}

### 1. Format patients' demographic (DX) data for JAGS run ---------

#number of patients
npat <- sapply(dx_list, nrow)

#number with known true state
npat_cancer_known <- sapply(dx_list, function(x){sum(!is.na(x$rp))})

#list of known true states
cancer_data_true <- matrix(0, nrow = J, ncol = max(npat_cancer_known))
for (i in 1:J){
  cancer_data_true[i,1:npat_cancer_known[i]] <- dx_list[[i]]$rp[!is.na(dx_list[[i]]$rp)]
}

#mean- and varianace- standardized age at diagnosis
for (i in 1:J){
  dx_list[[i]]$age_diag_std <- scale(dx_list[[i]]$age_diag)[,1]
  dx_list[[i]]$vol_std <- scale(dx_list[[i]]$vol)[,1]
}

# Covariate matrix for predictors of True gleasand grade
# Covariates can be any relevant predictors: here, age and vol
modmat_cancer <- array(0, dim = c(J, max(npat), 2))

for (i in 1:J){
  modmat_cancer[i,1:npat[i],1:2] <- as.matrix(cbind(dx_list[[i]]$age_diag_std, 
                                        dx_list[[i]]$vol_std))
}
  
npred_cancer <- rep(2,J)


### 2. Format PSA data for JAGS run ------------------
#number of observations
nobs_psa <- sapply(psa_list, nrow)

# observed PSA
log_psa_data <- psa_patient_index_map <- matrix(0, nrow = J, ncol = max(nobs_psa))
modmat_ranef_psa <- modmat_fixef_psa <- array(0, dim = c(J, max(nobs_psa), 2))

# covariate matrix for random effects in PSA mixed effects model
npred_ranef_psa <- npred_fixef_psa <- rep(2,J)

# covariate matrix for fixed effects in PSA mixed effects model

for (i in 1:J){
  age_tmp <- dx_list[[i]] %>% select(ID, age_diag_std)
  prosvol_tmp <- mri_list[[i]] %>% select(ID, mriprosvol) %>%
    group_by(ID) %>% summarise(prosvol_mean = mean(mriprosvol))
  prosvol_tmp$prosvol_std <- scale(prosvol_tmp$prosvol_mean)[,1]
  prosvol_tmp$prosvol_std[is.na(prosvol_tmp$prosvol_std)] <- 0
  prosvol_tmp <- prosvol_tmp %>% select(ID, prosvol_std)
  psa_list[[i]] <- psa_list[[i]] %>% left_join(age_tmp, by = "ID") %>%
    left_join(prosvol_tmp, by = "ID")
  psa_list[[i]]$prosvol_std[is.na(psa_list[[i]]$prosvol_std)] <- 0
}

for (i in 1:J){
  log_psa_data[i, 1:nobs_psa[i]] <- psa_list[[i]]$psa
  psa_patient_index_map[i, 1:nobs_psa[i]] <- psa_list[[i]]$ID
  modmat_ranef_psa[i,1:nobs_psa[i],1:2] <- as.matrix(cbind(rep(1,nobs_psa[i]), 
                                                           psa_list[[i]]$dxPSAdays))
  modmat_fixef_psa[i,1:nobs_psa[i],1:2] <- as.matrix(cbind(psa_list[[i]]$prosvol_std, psa_list[[i]]$age_diag_std))
}
data.check(condition=as.logical(sum(is.na(log_psa_data))==0), message="Missing PSA values. Email Yates; she will check orginal script.")
data.check(condition=as.logical(sum(is.na(modmat_ranef_psa))==0), message="Missing time since dx in PSA data. Email Yates; she will check orginal script.")

data.check(condition=as.logical(sum(is.na(modmat_fixef_psa))==0), message="Missing volumes in PSA data. Email Yates; she will check orginal script.")


# lmer fit to get initial value for covariance parameter
var_vec <- list(length = J)
for (i in 1:J){
  mod_lmer <-lmer(psa ~ prosvol_std + age_diag_std + (1 + dxPSAdays | ID), 
                  data = psa_list[[i]])
  var_vec[[i]] <- apply(coef(mod_lmer)$ID, 2, var)[1:npred_ranef_psa[i]]
  var_vec[[i]] <- var_vec[[i]][2:1]
}

### 5. Format biopsy result (reclassification) data for JAGS run -------------
#### IF NPC outcome was included in the model, would want to limit biopsy grade results to biopsies with cancer found
#make sure there are no missing biopsy grade outcomes

bxmri_list <- list(length = J)

for (i in 1:J){
  bx_tmp <- bx_list[[i]] %>% mutate(type = "bx")
  mri_tmp <- mri_list[[i]] %>% 
    mutate(type = "mri") %>% 
    select(ID, type, mripirads, dxMRIdays) %>% 
    rename(dxBXdays = dxMRIdays)
  bxmri_tmp <- full_join(bx_tmp, mri_tmp, by = c("ID", "type", "dxBXdays")) %>% 
    arrange(ID, dxBXdays)
  
  cur_ID <- 1
  cur_pirads <- NA
  for (k in 1:nrow(bxmri_tmp)){
    if (bxmri_tmp$type[k] == "mri"){
      cur_pirads <- bxmri_tmp$mripirads[k]
    } else if (bxmri_tmp$ID[k] == cur_ID){
      bxmri_tmp$mripirads[k] <- cur_pirads
    } else{
      cur_ID <- bxmri_tmp$ID[k]
      cur_pirads <- NA
    }
  }
  bxmri_list[[i]] <- bxmri_tmp %>% filter(type == "bx") %>%
    select(-c(type, bxpos))
}

gleason_to_pgg <- function(gleason) {
  if (sum(gleason) <= 6) {
    return(1)
  } else if (gleason[1] == 3){
    return (2)
  } else if (gleason[1] == 4 & gleason[2] == 3){
    return(3)
  } else return(4)
}

for (i in 1:J){
  bxmri_list[[i]]$pgg <- apply(bxmri_list[[i]][,2:3], 1, gleason_to_pgg)
}

# Define Cancer state argument
cancer_state <- matrix(0, nrow = J, ncol = max(npat - npat_cancer_known))
for (i in 1:J){
  cancer_state[i,1:((npat - npat_cancer_known)[i])] <- c(bxmri_list[[i]] %>% group_by(ID) %>%
    summarise(pgg_mean = round(mean(pgg))) %>%
    right_join(dx_list[[i]], by = "ID") %>%
    filter(is.na(rp)) %>% select(pgg_mean))$pgg_mean
}

#number of biopsies with GS outcome
n_pgg <- sapply(bxmri_list, nrow)

#GS outcomes for biopsies
PGG <- lapply(bxmri_list, function(x){x$pgg})
data.check(condition=as.logical(sum(is.na(PGG)) == 0), message="Missing biopsy results data. Email Yates; she will check orginal script.")
subj_pgg <- lapply(bxmri_list, function(x){x$ID})

npat_pgg <- n_pgg
pgg_data <- pgg_patient_index_map <- pgg_pirads_data <- pgg_pirads_data_m2 <- matrix(0, nrow = J, ncol = max(n_pgg))
for (i in 1:J){
  pgg_data[i,1:n_pgg[i]] <- PGG[[i]]
  pgg_patient_index_map[i,1:n_pgg[i]] <- subj_pgg[[i]]
}

# Covariates (beta) for Biopsy model: here, ratio of positive, natural spline of days until biopsy. 
modmat_pgg <- array(0, dim = c(J, max(n_pgg), 4))

for (i in 1:J){
  posratio_tmp <- bxmri_list[[i]]$corespos / bxmri_list[[i]]$corestaken
  posratio_tmp[is.na(posratio_tmp)] <- 0
  bxmri_list[[i]]$mripirads[is.na(bxmri_list[[i]]$mripirads)] <- 0
  
  pgg_pirads_data[i,1:n_pgg[i]] <- bxmri_list[[i]]$mripirads
  pgg_pirads_data_m2[i,1:n_pgg[i]] <- bxmri_list[[i]]$mripirads - 2
  
  #natural splines for continuous variables
  modmat_pgg[i,1:n_pgg[i], 1:4] <- as.matrix(cbind(posratio_tmp, ns(bxmri_list[[i]]$dxBXdays, 3)))
  npred_pgg <- rep(4,J)
}


### 6. Prior means for intercept in continuation ratio model ----------

logit <- function(x) log(x / (1-x))


# Latent cancer state
cancer_int1_mean <- cancer_int2_mean <- cancer_int3_mean <- rep(0,J)
for (i in 1:J){
  phat <- rep(0,4)
  for (k in 1:4){
    phat[k] <- sum(cancer_data_true[i,1:npat_cancer_known[i]] == k) / npat_cancer_known[i]
  }
  cancer_int1_mean[i] <- logit(phat[1])
  cancer_int2_mean[i] <- logit(phat[2] / (1 - phat[1]))
  cancer_int3_mean[i] <- logit(phat[3] / (1 - phat[1] - phat[2]))
}

# Biopsy 
pgg_int1_mean <- pgg_int2_mean <- pgg_int3_mean <- rep(0,J)
for (i in 1:J){
  phat_b <- rep(0,4)
  for (k in 1:4){
    phat_b[k] <- sum(PGG[[i]] == k) / length(PGG[[i]])
  }
  pgg_int1_mean[i] <- logit(phat_b[1])
  pgg_int2_mean[i] <- logit(phat_b[2] / (1 - phat_b[1]))
  pgg_int3_mean[i] <- logit(phat_b[3] / (1 - phat_b[1] - phat_b[2]))
}



