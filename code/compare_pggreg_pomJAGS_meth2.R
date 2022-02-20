## compare sequential model and proportional odds using patients under surgery (eta known)

rm(list=ls())
source('R-scripts/load_libs.R')
base.location <- "/Users/zitongwang/Downloads/prostate-active-surveillance-vDaan/" #"/users/rcoley/jhas-epic/psa/psa-new/" #"/Users/ryc/Documents/inhealth/prediction-model/automated/for-TIC/"
setwd(base.location)
location.of.data <- paste0(base.location, "data")
location.of.r.scripts <- paste0(base.location, "R-scripts")
location.of.generated.files <- paste0(base.location, "generated-files-etaknown-model/test")

name.of.pt.data <- "Demographic_6.15.csv" #"demographics with physician info.2015.csv"
name.of.bx.data <- "Biopsy_6.15.csv" #"Biopsy data_2015.csv"
name.of.psa.data <- "PSA_6.15.csv" #"PSA.2015.csv"
name.of.tx.data <- "Treatment_6.15.csv" #"treatment_2015.csv"

date.pull<-as.numeric(Sys.Date())
source(paste(location.of.r.scripts,"data-load-check-and-shaping.R",sep="/"))
#source(paste(location.of.r.scripts,"data-prep-for-jags-reform.R",sep="/"))
#source(paste(location.of.r.scripts,"argument-prep-for-jags-reform.R",sep="/"))
# source(paste(location.of.r.scripts,"JAGS-prediction-model-etaknown.R",sep="/"))

eta.info <- pt.data %>% dplyr::select(clinical_PTnum,subj, true.pgg)
psa.data.eta <- psa.data %>% left_join(eta.info)
psa.data.eta <- psa.data.eta %>% filter(!is.na(true.pgg))

# use only known eta patients
pt.data.eta1 <- pt.data %>% filter(!is.na(true.pgg))
pt.data.eta2 <- pt.data %>% filter(surg == 1)
pt.data.eta2$id[which(!(pt.data.eta2$id %in% pt.data.eta1$id))] # people went through surgery still has unknown pgg
pt.data.eta <- pt.data.eta1
npat <- dim(pt.data.eta)[1]
cancer_data <- pt.data.eta$true.pgg[!is.na(pt.data.eta$true.pgg)]
npat_cancer_known <- length(cancer_data)
#biospy result
pt.data.eta$dx.age.std <- scale(pt.data.eta$dx.age)
bx.full <- bx.full %>% left_join(eta.info)
bx.full.eta <- bx.full %>% filter(clinical_PTnum %in% pt.data.eta$clinical_PTnum)

bx.full.eta$pgg[bx.full.eta$npc==0 & !is.na(bx.full.eta$npc)]<-1
pgg.data0 <- bx.full.eta[bx.full.eta$bx.here==1 & !is.na(bx.full.eta$bx.here)
                        & !is.na(bx.full.eta$pgg)
                        &  bx.full.eta$time.int>0,]
pgg.data1 <- pgg.data0 %>% dplyr::select(clinical_PTnum, subj, bx.dt.num, ncs, mri, std.vol,
                                         true.pgg,pgg)
n_pgg1 <- dim(pgg.data1)[1]

pgg.data2 <- pgg.data0[!is.na(pgg.data0$bx.time.min),]
pgg.data2 <- pgg.data2 %>% dplyr::select(clinical_PTnum, subj, 
                                         bx.dt.num.min, ncs.min, mri.min, std.vol,
                                         true.pgg,pgg) %>% 
  dplyr::rename(bx.dt.num = bx.dt.num.min,
                ncs = ncs.min,
                mri = mri.min)
n_pgg2 <- dim(pgg.data2)[1]

pgg.data <- rbind(pgg.data1, pgg.data2)
#pgg.data$eta <- cancer_data[pgg.data$subj] # this is the same as pgg.data$true.pgg
# true pgg used for step functions in regression, pgg_data are measured pgg
pgg.data$eta_m2 <- ifelse(pgg.data$true.pgg-2 >=0 , 1, 0)
pgg.data$eta_m3 <- ifelse(pgg.data$true.pgg-3 >=0 , 1, 0)
pgg.data$eta_m4 <- ifelse(pgg.data$true.pgg-4 >=0 , 1, 0)

bx.date.knots <- attr(ns(pgg.data$bx.dt.num, 3), "knots")
bx.date.bknots <- attr(ns(pgg.data$bx.dt.num, 3), "Boundary.knots")
pgg.ncs.knots <- attr(ns(pgg.data$ncs, 2), "knots")
pgg.ncs.bknots <- attr(ns(pgg.data$ncs, 2), "Boundary.knots")

modmat_pgg <- as.matrix(cbind( ns(pgg.data$bx.dt.num, Boundary.knots = bx.date.bknots, knots = bx.date.knots),
                               ns(pgg.data$ncs , Boundary.knots = pgg.ncs.bknots, knots = pgg.ncs.knots),
                               pgg.data[, c("mri", "std.vol", "eta_m2", "eta_m3", "eta_m4")] ))
npred_pgg <- dim(modmat_pgg)[2]

colnames(modmat_pgg) <- c("ns(bx.dt.num,1)","ns(bx.dt.num,2)","ns(bx.dt.num,3)",
                     "ns(ncs,1)", "ns(ncs,2)", 
                     "mri", "std.vol", "eta_m2", "eta_m3", "eta_m4")
pgg_data <- pgg.data$pgg
npat_pgg <- n_pgg1 + n_pgg2
#The number of latent classes/ values of true cancer state
nlevel_cancer <- 4


### 1. Set up jags arguments
jags_data<-list(nlevel_cancer=nlevel_cancer, 
                pgg_data=pgg_data, npat_pgg=npat_pgg, 
                modmat_pgg=modmat_pgg, npred_pgg=npred_pgg
) 

### 2. Initialize model parameters

inits <- function() {
  pgg_intercept0 <- rnorm((nlevel_cancer-1),0,1)
  # pgg_slope <- rnorm((npred_pgg+(nlevel_cancer-1)),mean=0,sd=0.25)
  # pgg_slope[(npred_pgg+1):(npred_pgg+(nlevel_cancer-1))] <- 
  #   abs(pgg_slope[(npred_pgg+1):(npred_pgg+(nlevel_cancer-1))])
  pgg_slope <- rnorm((npred_pgg-(nlevel_cancer-1)),mean=0,sd=0.25)
  pgg_slope[(npred_pgg-(nlevel_cancer-1)+1):(npred_pgg)] <- 
    abs(pgg_slope[(npred_pgg+1):(npred_pgg+(nlevel_cancer-1))])
  
  list(
    pgg_intercept0=pgg_intercept0, pgg_slope=pgg_slope
  ) }


### 3. Define parameters to be tracked
params <- c("pgg_intercept", "pgg_slope")



### 4. Define other jags settings  
n.iter <- 10000; n.burnin <- 2500; n.thin <- 10; n.chains <- 1

### 5. Run
seed=2022
set.seed(seed)
outj<-jags(jags_data, inits=inits, parameters.to.save=params,
           model.file=paste("code/compare_pggreg_pomJAGS_meth2.txt"),
           n.thin=n.thin, n.chains=n.chains, n.burnin=n.burnin, n.iter=n.iter)
out<-outj$BUGSoutput
for(j in 1:length(out$sims.list)){
  write.csv(out$sims.list[[j]],
            paste(location.of.generated.files, "/jags-prediction-", names(out$sims.list)[j],"-",seed,".csv",sep=""))
  
}
get_stats <- function(x){a <- quantile(x, c(0.025, 0.975)); b<-mean(x); return(c(a[1], b, a[2]))}
jags_gamma_int <- 
  read_csv(paste0(location.of.generated.files,"/jags-prediction-pgg_intercept-2022.csv"))
apply(jags_gamma_int, 2, get_stats)
jags_gamma_int2 <- 
  read_csv(paste0("/Users/zitongwang/Downloads/prostate-active-surveillance-vDaan/generated-files-pgg-model","/jags-prediction-pgg_intercept-2022.csv"))
apply(jags_gamma_int2, 2, get_stats)

jags_gamma_slope <- 
  read_csv(paste0(location.of.generated.files,"/jags-prediction-pgg_slope-2022.csv"))
apply(jags_gamma_slope, 2, get_stats)

jags_gamma_slope2 <- 
  read_csv(paste0("/Users/zitongwang/Downloads/prostate-active-surveillance-vDaan/generated-files-etaknown-model","/jags-prediction-pgg_slope-2022.csv"))
apply(jags_gamma_slope2, 2, get_stats)


