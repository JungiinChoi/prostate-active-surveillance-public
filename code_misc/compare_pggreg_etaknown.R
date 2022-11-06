## compare sequential model and proportional odds using patients under surgery (eta known)

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
#source(paste(location.of.r.scripts,"data-prep-for-jags.R",sep="/"))
#source(paste(location.of.r.scripts,"argument-prep-for-jags.R",sep="/"))
# source(paste(location.of.r.scripts,"JAGS-prediction-model-etaknown.R",sep="/"))

eta.info <- pt.data %>% dplyr::select(clinical_PTnum,subj, true.pgg)
psa.data.eta <- psa.data %>% left_join(eta.info)
psa.data.eta <- psa.data.eta %>% filter(!is.na(true.pgg))

# use only known eta patients
pt.data.eta1 <- pt.data %>% filter(!is.na(true.pgg))
pt.data.eta2 <- pt.data %>% filter(surg == 1)
pt.data.eta2$id[which(!(pt.data.eta2$id %in% pt.data.eta1$id))] # people went through surgery still has unknown pgg
pt.data.eta <- pt.data.eta1

n <- dim(pt.data.eta)[1]
eta.data <- pt.data.eta$true.pgg[!is.na(pt.data.eta$true.pgg)]
n_eta_known <- length(eta.data)
#biospy result
pt.data.eta$dx.age.std <- scale(pt.data.eta$dx.age)
bx.full <- bx.full %>% left_join(eta.info)
bx.full.eta <- bx.full %>% filter(clinical_PTnum %in% pt.data.eta$clinical_PTnum)
bx.full.eta$pgg[bx.full.eta$npc==0 & !is.na(bx.full.eta$npc)]<-1
pgg.data <- bx.full.eta[bx.full.eta$bx.here==1 & !is.na(bx.full.eta$bx.here)
                        & !is.na(bx.full.eta$pgg)
                        &  bx.full.eta$time.int>0,]
n_pgg <- dim(pgg.data)[1]
PGG <- as.numeric(pgg.data$pgg)
subj_pgg <- pgg.data$subj
V.PGG.ID <- pgg.data$clinical_PTnum
V.PGG.data <- as.matrix(cbind(pgg.data$bx.dt.num, pgg.data$ncs, pgg.data$mri, pgg.data$std.vol))
d.V.PGG <- dim(V.PGG.data)[2]
pgg.data2 <- pgg.data[!is.na(pgg.data$bx.time.min),]
n_pgg2 <- dim(pgg.data2)[1]
PGG2 <- rep(1,n_pgg2)
subj_pgg2 <- pgg.data2$subj
V.PGG2.ID <- pgg.data2$clinical_PTnum
V.PGG2.data <- as.matrix(cbind(pgg.data2$bx.dt.num.min, pgg.data2$ncs.min,
                               pgg.data2$mri.min, pgg.data2$std.vol))

#combine all outcomes
n_pgg <- n_pgg + n_pgg2
PGG <- c(PGG, PGG2)
subj_pgg <- c(subj_pgg, subj_pgg2)
data.check(condition=as.logical(sum(is.na(PGG))==0), message="Missing biopsy results data. Email Yates; she will check orginal script.")
V.PGG.ID <- c(V.PGG.ID, V.PGG2.ID)
V.PGG.data <- rbind(V.PGG.data, V.PGG2.data)
bx.date.knots <- attr(ns(V.PGG.data[,1], 3), "knots")
bx.date.bknots <- attr(ns(V.PGG.data[,1], 3), "Boundary.knots")
pgg.ncs.knots <- attr(ns(V.PGG.data[,2], 2), "knots")
pgg.ncs.bknots <- attr(ns(V.PGG.data[,2], 2), "Boundary.knots")

V.PGG.data <- as.matrix(cbind( ns(V.PGG.data[,1], 3),
                               ns(V.PGG.data[,2], 2),
                               V.PGG.data[,3:d.V.PGG] ))
d.V.PGG <- dim(V.PGG.data)[2]
library(MASS)
pom.PGG.data <- as.data.frame(V.PGG.data)
pgg.data.all <- rbind(pgg.data, pgg.data2)
pom.PGG.data$eta_m2 <- ifelse(pgg.data.all$true.pgg-2 >=0 , 1, 0)
pom.PGG.data$eta_m3 <- ifelse(pgg.data.all$true.pgg-3 >=0 , 1, 0)
pom.PGG.data$eta_m4 <- ifelse(pgg.data.all$true.pgg-4 >=0 , 1, 0)
pom.PGG.data <- as.matrix(pom.PGG.data)
pom_fit <- polr(factor(PGG) ~ pom.PGG.data)

glmdata <- as.data.frame(V.PGG.data)
colnames(glmdata) <- c("ns(bx.dt.num,1)","ns(bx.dt.num,2)","ns(bx.dt.num,3)",
                       "ns(ncs,1)", "ns(ncs,2)", 
                       "mri", "std.vol")
glmdata$clinical_PTnum <- V.PGG.ID
glmdata$eta_m2 <- ifelse(pgg.data.all$true.pgg-2 >=0 , 1, 0)
glmdata$eta_m3 <- ifelse(pgg.data.all$true.pgg-3 >=0 , 1, 0)
glmdata$eta_m4 <- ifelse(pgg.data.all$true.pgg-4 >=0 , 1, 0)
glmdata$PGG <- PGG

glm_fit1<-glm(PGG_lev1 ~., family = binomial, 
          data = glmdata %>% 
            mutate(PGG_lev1 = factor(ifelse(PGG == 1, 1, 0), levels=c(0,1))) %>% 
            dplyr::select(-clinical_PTnum, -PGG))
glm_fit2<-glm(PGG_lev2 ~., family = binomial, 
          data = glmdata %>% 
            filter(PGG > 1) %>% 
            mutate(PGG_lev2 = factor(ifelse(PGG == 2, 1, 0), levels = c(0,1))) %>% 
            dplyr::select(-clinical_PTnum, -PGG))
glm_fit3<-glm(PGG_lev3 ~., family = binomial, 
          data = glmdata %>% 
            filter(PGG > 2) %>% 
            mutate(PGG_lev3 = factor(ifelse(PGG == 3, 1, 0), levels=c(0,1))) %>% 
            dplyr::select(-clinical_PTnum, -PGG))
p_y1_given_x0 <- predict(glm_fit1, type = "response")
p_y2_given_x0ge1 <- predict(glm_fit2, newdata = glmdata, type = "response")
p_y3_given_x0ge2 <- predict(glm_fit3, newdata = glmdata, type = "response")

glmdata$p1 <- p_y1_given_x0
glmdata$p2<- (1-p_y1_given_x0) * p_y2_given_x0ge1
glmdata$p3 <- (1-p_y1_given_x0)*(1-p_y2_given_x0ge1) * p_y3_given_x0ge2
glmdata$p4 <- 1-glmdata$p1-glmdata$p2-glmdata$p3
table(glmdata$p1+glmdata$p2 + glmdata$p3 + glmdata$p4)

#overall comparison
obs1 <- sum(glmdata$PGG == 1); exp1 <- sum(glmdata$p1)
obs2 <- sum(glmdata$PGG == 2); exp2 <- sum(glmdata$p2)
obs3 <- sum(glmdata$PGG == 3); exp3 <- sum(glmdata$p3)
obs4 <- sum(glmdata$PGG == 4); exp4 <- sum(glmdata$p4)
(obs1-exp1)^2/exp1 + (obs2-exp2)^2/exp2 + (obs3-exp3)^2/exp3 + (obs4-exp4)^2/exp4
##pgg=3 overestimate and pgg=4underestimate

pom_preds <- predict(pom_fit, type = "probs")
glmdata$pom_p1 <- pom_preds[,1]; glmdata$pom_p2 <- pom_preds[,2]
glmdata$pom_p3 <- pom_preds[,3]; glmdata$pom_p4 <- pom_preds[,4]
pom_exp1 <- sum(glmdata$pom_p1); pom_exp2 <- sum(glmdata$pom_p2)
pom_exp3 <- sum(glmdata$pom_p3); pom_exp4 <- sum(glmdata$pom_p4)
(obs1-pom_exp1)^2/pom_exp1 + (obs2-pom_exp2)^2/pom_exp2 + 
  (obs3-pom_exp3)^2/pom_exp3 + (obs4-pom_exp4)^2/pom_exp4

#by subgroup mri
## stratefy by MRI
table(glmdata$mri)
table(pom.PGG.data$mri)
get_obs_exp <- function(data){
  obs1 <- sum(data$PGG == 1); exp1 <- sum(data$p1)
  obs2 <- sum(data$PGG == 2); exp2 <- sum(data$p2)
  obs3 <- sum(data$PGG == 3); exp3 <- sum(data$p3)
  obs4 <- sum(data$PGG == 4); exp4 <- sum(data$p4)
  chisq_glm <- (obs1-exp1)^2/exp1 + (obs2-exp2)^2/exp2 + (obs3-exp3)^2/exp3 + (obs4-exp4)^2/exp4
  pom_exp1 <- sum(data$pom_p1); pom_exp2 <- sum(data$pom_p2)
  pom_exp3 <- sum(data$pom_p3); pom_exp4 <- sum(data$pom_p4)
  chisq_pom <- (obs1-pom_exp1)^2/pom_exp1 + (obs2-pom_exp2)^2/pom_exp2 + 
    (obs3-pom_exp3)^2/pom_exp3 + (obs4-pom_exp4)^2/pom_exp4
  return(data.frame(obs = c(obs1, obs2, obs3, obs4),
                    exp_glm = c(exp1, exp2, exp3, exp4),
                    exp_pom = c(pom_exp1,pom_exp2,pom_exp3,pom_exp4),
                    chisq_glm = rep(chisq_glm, 4),
                    chisq_pom = rep(chisq_pom,4)))
}
get_obs_exp(data = glmdata)
get_obs_exp(data = glmdata[glmdata$mri == 1,])
get_obs_exp(data = glmdata[glmdata$mri == 0,])

