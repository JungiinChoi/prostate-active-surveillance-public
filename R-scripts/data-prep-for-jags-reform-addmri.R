#reformat data prep code -----
#K -> nlevel_cancer
#K.bin  -> nlevel_cancer_bin
#n -> npat
#eta.data -> cancer_data
#n_eta_known -> npat_cancer_known
#V.eta.data -> modmat_cancer
#d.V.ETA -> npred_cancer
#n_obs_psa -> nobs_psa
#Y -> log_psa_data
#subj_psa -> psa_patient_index_map
#Z.data -> modmat_ranef_psa
#X.data -> modmat_fixef_psa
#d.Z -> npred_ranef_psa
#d.X -> npred_fixef_psa
#PGG -> pgg_data
#n_pgg -> npat_pgg
#subj_pgg -> pgg_patient_index_map
#V.PGG.data -> modmat_pgg
#d.V.PGG -> npred_pgg

#------------------------------
### 1. Format pt-level data for JAGS run
#number of patients
npat <- dim(pt.data)[1]

#list of known true states
cancer_data <- pt.data$true.pgg[!is.na(pt.data$true.pgg)]

#number with known true state
npat_cancer_known <- length(cancer_data)

#mean- and varianace- standardized age at diagnosis
pt.data$dx.age.std <- scale(pt.data$dx.age)

#covariate matrix for predictors of PGG
modmat_cancer <- as.matrix(cbind(pt.data$dx.age.std, pt.data$lr.vol))
npred_cancer<-dim(modmat_cancer)[2]



### 2. Format PSA data for JAGS run
#number of observations
nobs_psa <- dim(psa.data)[1]

#observed PSA
log_psa_data <- psa.data$log.psa
data.check(condition=as.logical(sum(is.na(log_psa_data))==0), message="Missing PSA values. Email Yates; she will check orginal script.")

#list of patients
psa_patient_index_map <- psa.data$subj

#covariate matrix for random effects
modmat_ranef_psa <- as.matrix(cbind(rep(1,nobs_psa), psa.data$time.since.dx))
npred_ranef_psa <- dim(modmat_ranef_psa)[2]
data.check(condition=as.logical(sum(is.na(modmat_ranef_psa))==0), message="Missing time since dx in PSA data. Email Yates; she will check orginal script.")

#covariate matrix for fixed effects
#take dx age from pt.data and put into psa.data
psa.data$dx.age.std<-vector(length=nobs_psa)
for(i in 1:npat){
  psa.data$dx.age.std[psa.data$subj==pt.data$subj[i]]<-pt.data$dx.age.std[i]}

#define covariate matrix
modmat_fixef_psa <- as.matrix(cbind(psa.data$std.vol, scale(psa.data$dx.age.std)))
npred_fixef_psa <- dim(modmat_fixef_psa)[2]
data.check(condition=as.logical(sum(is.na(modmat_fixef_psa))==0), message="Missing volumes in PSA data. Email Yates; she will check orginal script.")

#lmer fit to get starting value for covariance parameter
mod_lmer<-lmer(log.psa~ std.vol + dx.age.std + (1+ time.since.dx |subj), data=psa.data)
var_vec <- apply(coef(mod_lmer)$subj, 2, var)[1:npred_ranef_psa]
names(var_vec)
#make sure order of variance estimates is correct
index.intercept<-c(1:2)[names(var_vec)=="(Intercept)"]
index.time<-c(1:2)[names(var_vec)=="time.since.dx"]
var_vec <- c(var_vec[index.intercept], var_vec[index.time])


### 5. Format biopsy result (reclassification) data for JAGS run
#### IF NPC outcome was included in the model, would want to limir biopsy grade results to biopsies with cancer found
#make sure there are no missing biopsy grade outcomes
bx.full$pgg[bx.full$npc==0 & !is.na(bx.full$npc)]<-1

#subset data to intervals with biopsy occurring
pgg.data <- bx.full[bx.full$bx.here==1 & !is.na(bx.full$bx.here)
                    & !is.na(bx.full$pgg)
                    &  bx.full$time.int>0,]

#number of biopsies with GS outcome
n_pgg <- dim(pgg.data)[1]

#GS outcomes for biopsies
PGG <- as.numeric(pgg.data$pgg)

#unique subj identifier
subj_pgg <- pgg.data$subj

#covariate matrix
V.PGG.data <- as.matrix(cbind(pgg.data$bx.dt.num, pgg.data$ncs, pgg.data$mri, pgg.data$std.vol))
d.V.PGG <- dim(V.PGG.data)[2]

#for intervals with second biopsies
pgg.data2 <- pgg.data[!is.na(pgg.data$bx.time.min),]
                      #& !is.na(pgg.data$npc.min) & pgg.data$npc.min>=0,]
n_pgg2 <- dim(pgg.data2)[1]
PGG2 <- rep(1,n_pgg2)
subj_pgg2 <- pgg.data2$subj
V.PGG2.data <- as.matrix(cbind(pgg.data2$bx.dt.num.min, pgg.data2$ncs.min,
                               pgg.data2$mri.min, pgg.data2$std.vol))

#combine all outcomes
npat_pgg <- n_pgg + n_pgg2
pgg_data <- c(PGG, PGG2)
pgg_patient_index_map <- c(subj_pgg, subj_pgg2)
data.check(condition=as.logical(sum(is.na(pgg_data))==0), message="Missing biopsy results data. Email Yates; she will check orginal script.")

V.PGG.data <- rbind(V.PGG.data, V.PGG2.data)
#natural splines for continuous variables
bx.date.knots <- attr(ns(V.PGG.data[,1], 3), "knots")
bx.date.bknots <- attr(ns(V.PGG.data[,1], 3), "Boundary.knots")
pgg.ncs.knots <- attr(ns(V.PGG.data[,2], 2), "knots")
pgg.ncs.bknots <- attr(ns(V.PGG.data[,2], 2), "Boundary.knots")

modmat_pgg <- as.matrix(cbind( ns(V.PGG.data[,1], 3),
                               ns(V.PGG.data[,2], 2),
                               V.PGG.data[,3:d.V.PGG] ))
npred_pgg <- dim(modmat_pgg)[2]
data.check(condition=as.logical(sum(is.na(V.PGG.data))==0), message="Missing biopsy data. Email Yates; she will check orginal script.")

## MRI data
mri_data_615 <- read_csv("data/MRI_6.15.csv")
mri_data <- mri_data_615 %>% dplyr::select(clinical_ptnum, pirads) %>% group_by(clinical_ptnum) %>% summarize(pirads_max = max(pirads, na.rm = T))

mri_data$pirads_max <- ifelse(mri_data$pirads_max == "-Inf", NA, mri_data$pirads_max)
mri_data$id<-mri_data$clinical_ptnum
mri_data$subj<-rep(0,dim(mri_data)[1])
for(i in 1:dim(mri_data)[1]){
  mri_data$subj[mri_data$id==pt.data$id[i]] <- i}
#mri_data$id[mri_data$subj == 0]
#sum(mri_data$id[mri_data$subj == 0] %in% unique(pt.data$clinical_PTnum))
mri_data <- mri_data[mri_data$subj != 0,]
n_mri <- dim(mri_data)[1]

mri_data_complete <- mri_data[complete.cases(mri_data$pirads_max),]
pirads_data <- mri_data_complete$pirads_max
pirads_patient_index_map <- mri_data_complete$subj
npat_pirads <- dim(mri_data_complete)[1]
npred_pirads<-0

nobs_pirads <- dim(mri_data_complete)[1]
