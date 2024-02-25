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
J <- 3

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

#list of known true states
cancer_data_true <- sapply(dx_list, function(x){x$rp[!is.na(x$rp)]})


dx_list[[3]]$rp == "NA"
for (i in 1:J){
  cancer_data[[i]] <- (dx_list[[i]]$rp)[order(is.na(dx_list[[i]]$rp))]
}

#number with known true state
npat_cancer_known <- sapply(dx_list, function(x){sum(!is.na(x$rp))})

#mean- and varianace- standardized age at diagnosis
for (i in 1:J){
  dx_list[[i]]$age_diag_std <- scale(dx_list[[i]]$age_diag)[,1]
  dx_list[[i]]$vol_std <- scale(dx_list[[i]]$vol)[,1]
}

# Covariate matrix for predictors of True gleasand grade
# Covariates can be any relevant predictors: here, age and vol
modmat_cancer <- list(length = J)

for (i in 1:J){
  modmat_cancer[[i]] <- as.matrix(cbind(dx_list[[i]]$age_diag_std, 
                                        dx_list[[i]]$vol_std))
}
  
npred_cancer <- sapply(modmat_cancer, ncol)


### 2. Format PSA data for JAGS run ------------------
#number of observations
nobs_psa <- sapply(psa_list, nrow)

#observed PSA
log_psa_data <- lapply(psa_list, function(x){x$psa})
data.check(condition=as.logical(sum(is.na(log_psa_data))==0), message="Missing PSA values. Email Yates; she will check orginal script.")

#list of patients
psa_patient_index_map <- lapply(psa_list, function(x){x$ID})

# covariate matrix for random effects in PSA mixed effects model
modmat_ranef_psa <- list(length = J)
for (i in 1:J){
  modmat_ranef_psa[[i]] <- as.matrix(cbind(rep(1,nobs_psa[i]), 
                                           psa_list[[i]]$dxPSAdays))
}
npred_ranef_psa <- sapply(modmat_ranef_psa, ncol)
data.check(condition=as.logical(sum(is.na(modmat_ranef_psa))==0), message="Missing time since dx in PSA data. Email Yates; she will check orginal script.")

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

modmat_fixef_psa <- lapply(psa_list, function(x){as.matrix(cbind(x$prosvol_std, x$age_diag_std))})
npred_fixef_psa <- sapply(modmat_fixef_psa, ncol)
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
cancer_state = list(length = J)
for (i in 1:J){
  cancer_state[[i]] <- bxmri_list[[2]] %>% group_by(ID) %>%
    summarise(pgg_mean = round(mean(pgg))) %>%
    right_join(dx_list[[i]], by = "ID") %>%
    filter(is.na(rp)) %>%
    select(pgg_mean)
}

#number of biopsies with GS outcome
n_pgg <- sapply(bxmri_list, nrow)

#GS outcomes for biopsies
PGG <- lapply(bxmri_list, function(x){x$pgg})
data.check(condition=as.logical(sum(is.na(PGG)) == 0), message="Missing biopsy results data. Email Yates; she will check orginal script.")
subj_pgg <- lapply(bxmri_list, function(x){x$ID})

npat_pgg <- n_pgg
pgg_data <- PGG
pgg_patient_index_map <- subj_pgg

# Covariates (beta) for Biopsy model: here, ratio of positive, natural spline of days until biopsy. 

pgg_pirads_data <- list(length = J)
pgg_pirads_data_m2 <- list(length = J)
modmat_pgg <- list(length = J)
for (i in 1:J){
  posratio_tmp <- bxmri_list[[i]]$corespos / bxmri_list[[i]]$corestaken
  posratio_tmp[is.na(posratio_tmp)] <- 0
  bxmri_list[[i]]$mripirads[is.na(bxmri_list[[i]]$mripirads)] <- 0
  pgg_pirads_data[[i]] <- bxmri_list[[i]]$mripirads
  
  pgg_pirads_data_m2[[i]] <- pgg_pirads_data[[i]] - 2
  
  #natural splines for continuous variables
  modmat_pgg[[i]] <- as.matrix(cbind(posratio_tmp, ns(bxmri_list[[i]]$dxBXdays, 3)))
  npred_pgg <- sapply(modmat_pgg, ncol)
}

data.check(condition=as.logical(sum(is.na(V.PGG.data))==0), message="Missing biopsy data. Email Yates; she will check orginal script.")


### 6. Prior means for intercept in continuation ratio model ----------

logit <- function(x) log(x / (1-x))


# Latent cancer state
cancer_int1_mean <- cancer_int2_mean <- cancer_int3_mean <- rep(0,J)
for (i in 1:J){
  phat <- rep(0,4)
  for (k in 1:4){
    phat[k] <- sum(cancer_data_true[[i]] == k) / length(cancer_data_true[[i]])
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



