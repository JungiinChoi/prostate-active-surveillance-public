########### SEED/GROUP TO MASK IS EITHER SET HERE OR IN EXAMPLE-CV-TO-RUN.R
#### This should be adjusted by the user depending on how seed, masking set 
# SEED <- as.numeric(Sys.getenv("SGE_TASK_ID"))
# # For GAP3 analysis, make this seed=1,...,K and define grouping variable with values 1,...,K
# to.mask<- SEED
###########



# Yates Coley
# rycoley@gmail.com
# 2019 January 08, modified 2019 January 30
# This script will take tidied/shaped AS data and perform additional manipulation for the JAGS model
# I've put a lot of data checks in here as fail-safes, but there shouldn't be any missing values after running the previous script.
# This updated from data-prep-for-jags.R to accomodated K-fold CV 

#### WORKFLOW
### 1. Format pt-level data for JAGS run
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

# load data files

dx_data <- read.csv(paste0(location.of.data,"/dx.csv"))
bx_data <- read.csv(paste0(location.of.data,"/bx.csv"))
psa_data <- read.csv(paste0(location.of.data,"/psa.csv"))
mri_data <- read.csv(paste0(location.of.data,"/mri.csv"))

# ID preprocess
id_mapping <- function(x){which(dx_data$id == x)}
dx_data$ID <- 1:nrow(dx_data)
bx_data$ID <- sapply(bx_data$id, id_mapping)
psa_data$ID <- sapply(psa_data$id, id_mapping)
mri_data$ID <- sapply(mri_data$id, id_mapping)

### 1. Format patients' demographic (DX) data for JAGS run ---------

#number of patients
npat <- nrow(dx_data)

#list of known true states
cancer_data_true <- dx_data$rp[!is.na(dx_data$rp)]

cancer_data <- (dx_data$rp)[order(is.na(dx_data$rp))]

#number with known true state
npat_cancer_known <- sum(!is.na(dx_data$rp))

#mean- and varianace- standardized age at diagnosis

dx_data$age_diag_std <- scale(dx_data$age_diag)[,1]
dx_data$vol_std <- scale(dx_data$vol)[,1]


# Covariate matrix for predictors of True gleasand grade
# Covariates can be any relevant predictors: here, age and vol

dx_data$vol_std[is.na(dx_data$vol_std)] <- mean(dx_data$vol_std, na.rm = TRUE)
dx_data$age_diag_std[is.na(dx_data$age_diag_std)] <- mean(dx_data$age_diag_std, na.rm = TRUE)
modmat_cancer <- as.matrix(cbind(dx_data$age_diag_std, 
                                 dx_data$vol_std))

npred_cancer <- ncol(modmat_cancer)


### 2. Format PSA data for JAGS run ------------------
#number of observations
nobs_psa <- nrow(psa_data)

#observed PSA
log_psa_data <- psa_data$psa
data.check(condition=as.logical(sum(is.na(log_psa_data))==0), message="Missing PSA values. Email Yates; she will check orginal script.")

#list of patients
psa_patient_index_map <- psa_data$ID

# covariate matrix for random effects in PSA mixed effects model
modmat_ranef_psa <- as.matrix(cbind(rep(1,nobs_psa), 
                                         psa_data$dxPSAdays))

npred_ranef_psa <- ncol(modmat_ranef_psa)
data.check(condition=as.logical(sum(is.na(modmat_ranef_psa))==0), message="Missing time since dx in PSA data. Email Yates; she will check orginal script.")

# covariate matrix for fixed effects in PSA mixed effects model


age_tmp <- dx_data %>% select(ID, age_diag_std)
prosvol_tmp <- mri_data %>% select(ID, mriprosvol) %>%
  group_by(ID) %>% summarise(prosvol_mean = mean(mriprosvol))
prosvol_tmp$prosvol_std <- scale(prosvol_tmp$prosvol_mean)[,1]
prosvol_tmp$prosvol_std[is.na(prosvol_tmp$prosvol_std)] <- 0
prosvol_tmp <- prosvol_tmp %>% select(ID, prosvol_std)
psa_data <- psa_data %>% left_join(age_tmp, by = "ID") %>%
  left_join(prosvol_tmp, by = "ID")
psa_data$prosvol_std[is.na(psa_data$prosvol_std)] <- 0

psa_data$age_diag_std[is.na(psa_data$age_diag_std)] <- 0
modmat_fixef_psa <- as.matrix(cbind(psa_data$prosvol_std, psa_data$age_diag_std))
npred_fixef_psa <- ncol(modmat_fixef_psa)
data.check(condition=as.logical(sum(is.na(modmat_fixef_psa))==0), message="Missing volumes in PSA data. Email Yates; she will check orginal script.")


# lmer fit to get initial value for covariance parameter
psa_data_i <- as.data.frame(scale(psa_data))
mod_lmer <-lmer(psa ~ prosvol_std + age_diag_std + (1 + dxPSAdays | ID), 
                data = psa_data_i)
var_vec <- apply(coef(mod_lmer)$ID, 2, var)[1:npred_ranef_psa]
var_vec <- var_vec[2:1]


### 5. Format biopsy result (reclassification) data for JAGS run -------------
#### IF NPC outcome was included in the model, would want to limit biopsy grade results to biopsies with cancer found
#make sure there are no missing biopsy grade outcomes
bx_tmp <- bx_data %>% mutate(type = "bx")
mri_tmp <- mri_data %>% 
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
bxmri_data <- bxmri_tmp %>% filter(type == "bx") %>%
  select(-c(type, bxpos))

gleason_to_pgg <- function(gleason) {
  if (sum(gleason) <= 6) {
    return(1)
  } else if (gleason[1] == 3){
    return (2)
  } else if (gleason[1] == 4 & gleason[2] == 3){
    return(3)
  } else return(4)
}

bxmri_data$pgg <- apply(bxmri_data[,2:3], 1, gleason_to_pgg)


# Define Cancer state argument

  cancer_state <- bxmri_data %>% group_by(ID) %>%
    summarise(pgg_mean = round(mean(pgg))) %>%
    right_join(dx_data, by = "ID") %>%
    filter(is.na(rp)) %>%
    select(pgg_mean)

  
#number of biopsies with GS outcome
n_pgg <- nrow(bxmri_data)

#GS outcomes for biopsies
PGG <- bxmri_data$pgg
data.check(condition=as.logical(sum(is.na(PGG)) == 0), message="Missing biopsy results data. Email Yates; she will check orginal script.")
subj_pgg <- bxmri_data$ID

npat_pgg <- n_pgg
pgg_data <- PGG
pgg_patient_index_map <- subj_pgg

# Covariates (beta) for Biopsy model: here, ratio of positive, natural spline of days until biopsy. 

  posratio_tmp <- bxmri_data$corespos / bxmri_data$corestaken
  posratio_tmp[is.na(posratio_tmp)] <- 0
  bxmri_data$mripirads[is.na(bxmri_data$mripirads)] <- 0
  pgg_pirads_data <- bxmri_data$mripirads
  
  pgg_pirads_data_m2 <- pgg_pirads_data - 2
  
  #natural splines for continuous variables
  modmat_pgg <- as.matrix(cbind(posratio_tmp, ns(bxmri_data$dxBXdays, df = 3)))
  npred_pgg <- ncol(modmat_pgg)

### 6. Format MRI data for JAGS run -------------
## formate MRI data for outcome model
  

pirads_data <- ifelse(mri_data$mripirads %in% c(1, 2), 1,
                      ifelse(mri_data$mripirads %in% c(3), 2, 3))  
  
pirads_patient_index_map <- mri_data$ID
npat_pirads <- nrow(mri_data)
  

### 7. Prior means for intercept in continuation ratio model ----------

logit <- function(x) log(x / (1-x))


# Latent cancer state

  phat <- rep(0,4)
  for (k in 1:4){
    phat[k] <- sum(cancer_data_true == k) / length(cancer_data_true)
  }
  cancer_int1_mean <- logit(phat[1])
  cancer_int2_mean <- logit(phat[2] / (1 - phat[1]))
  cancer_int3_mean <- logit(phat[3] / (1 - phat[1] - phat[2]))


# Biopsy 

  phat_b <- rep(0,4)
  for (k in 1:4){
    phat_b[k] <- sum(PGG == k) / length(PGG)
  }
  pgg_int1_mean <- logit(phat_b[1])
  pgg_int2_mean <- logit(phat_b[2] / (1 - phat_b[1]))
  pgg_int3_mean <- logit(phat_b[3] / (1 - phat_b[1] - phat_b[2]))

  
# Pirads
  pirads_int1_mean <- pirads_int2_mean <- 0 
  
  phat_p <- rep(0,3)
  for (k in 1:3){
    phat_p[k] <- sum(pirads_data == k) / length(pirads_data)
  }
  pirads_int1_mean <- logit(phat_p[1])
  pirads_int2_mean <- logit(phat_p[2] / (1 - phat_p[1]))
  
  