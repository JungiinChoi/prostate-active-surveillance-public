#Yates Coley
#rycoley@gmail.com
#Last updated August 24, 2017
#Annotations updated December 21, 2017
#This script will take tidied/shaped AS data and perform additional shaping for JAGS
#I've put a lot of data checks in here as fail-safes, but there shouldn't be any missing values after running the previous script.


#### WORKFLOW
### 1. Format pt-level data for JAGS run
### 2. Format PSA data for JAGS run
### 3. Format biopsy received data for JAGS run #OMITTED
### 4. Format NPC data for JAGS run #OMITTED
### 5. Format biopsy result (reclassification) data for JAGS run
### 6. Format MPC data for JAGS run #OMITTED
### 7. Format surgery received data for JAGS run #OMITTED
### 8. Format cumulative laterality data for JAGS run #OMITTED

#### OMITTED sections: Including biopsy and treatment/surgery decisions and additional biopsy findings (NPC, MPC, LAT) did not improve model estimation, so they were not included in the final model.
#Code for these outcomes is commented out here but provided so that the user could see how variables would be incorporated


### 1. Format pt-level data for JAGS run
#number of patients
n <- dim(pt.data)[1]

#list of known true states
eta.data <- pt.data$true.pgg[!is.na(pt.data$true.pgg)]

#number with known true state
n_eta_known <- length(eta.data)

#mean- and varianace- standardized age at diagnosis
pt.data$dx.age.std <- scale(pt.data$dx.age)

#covariate matrix for predictors of PGG
V.ETA.data <- as.matrix(cbind(pt.data$dx.age.std, pt.data$lr.vol))
d.V.ETA<-dim(V.ETA.data)[2]



### 2. Format PSA data for JAGS run
#number of observations
n_obs_psa <- dim(psa.data)[1]

#observed PSA
Y <- psa.data$log.psa
data.check(condition=as.logical(sum(is.na(Y))==0), message="Missing PSA values. Email Yates; she will check orginal script.")

#list of patients
subj_psa <- psa.data$subj

#covariate matrix for random effects
Z.data <- as.matrix(cbind(rep(1,n_obs_psa), psa.data$time.since.dx))
d.Z <- dim(Z.data)[2]
data.check(condition=as.logical(sum(is.na(Z.data))==0), message="Missing time since dx in PSA data. Email Yates; she will check orginal script.")

#covariate matrix for fixed effects
#take dx age from pt.data and put into psa.data
psa.data$dx.age.std<-vector(length=n_obs_psa)
for(i in 1:n){
  psa.data$dx.age.std[psa.data$subj==pt.data$subj[i]]<-pt.data$dx.age.std[i]}

#define covariate matrix
X.data <- as.matrix(cbind(psa.data$std.vol, scale(psa.data$dx.age.std)))
d.X <- dim(X.data)[2]
data.check(condition=as.logical(sum(is.na(X.data))==0), message="Missing volumes in PSA data. Email Yates; she will check orginal script.")

#lmer fit to get starting value for covariance parameter
mod_lmer<-lmer(log.psa~ std.vol + dx.age.std + (1+ time.since.dx |subj), data=psa.data)
var_vec <- apply(coef(mod_lmer)$subj, 2, var)[1:d.Z]
names(var_vec)
#make sure order of variance estimates is correct
index.intercept<-c(1:2)[names(var_vec)=="(Intercept)"]
index.time<-c(1:2)[names(var_vec)=="time.since.dx"]
var_vec <- c(var_vec[index.intercept], var_vec[index.time])



### 3. Format biopsy received data for JAGS run

### BX HERE #Decision to get biopsy in an annual interval
##subset data to post-dx
#bx.data <- bx.full[bx.full$time.int>0 & !is.na(bx.full$bx.here),]

##number of biopsies
#n_bx <- dim(bx.data)[1]

##binary outcome, presence of biopsy
#BX <- as.numeric(bx.data$bx.here)

##unique pt id for each bx
#subj_bx <- bx.data$subj
#data.check(condition=as.logical(sum(is.na(BX))==0), message="Missing biopsy data. Email Yates; she will check orginal script.")

##assume patients with missing data had 0 prior biopsies without cancer found
#bx.data$prop.npc0.start[is.na(bx.data$prop.npc0.start)]<-0

##covariate matrix
#U.BX.data <- as.matrix( cbind( rep(1, n_bx),
#                               ns(bx.data$time.int, 3),
#                               ns(bx.data$int.dt.num, 3),
#                               ns(bx.data$int.age, 3),
#                               ns(bx.data$num.prev.bx.start, 3),
#                               ns(bx.data$max.prev.mpc.start, 2),
#                               bx.data$psa.pres.bx,
#                               bx.data$std.vol,
#                               ns(bx.data$prop.npc0.start,2),
#                               ns(bx.data$freq.bx.start,2)
#                               ))
#d.U.BX<-dim(U.BX.data)[2]
#data.check(condition=as.logical(sum(is.na(U.BX.data))==0), message="Missing values in biopsy data. Email Yates; she will check orginal script.")


## Regression fit below was examined to select covariates
##### In future, use lasso approach? how should I select flexibility for continuous variables?
#mod<-glm(BX ~ ns(bx.data$time.int, 3) +
 #        ns(bx.data$int.dt.num, 3) +
#         ns(bx.data$int.age, 3) +
#         ns(bx.data$num.prev.bx.start, 3) +
#         ns(bx.data$max.prev.mpc.start, 2) +
#         #ns(bx.data$max.prev.npc.start, 3) +
#         bx.data$psa.pres.bx +# bx.data$psa.traj.bx +
#         bx.data$std.vol +
#         #bx.data$prev.vol.rc.start +
#         ns(bx.data$prop.npc0.start,2) +
#         ns(bx.data$freq.bx.start,2) +
#           bx.data$prev.lat.start,
#         family="binomial")
#summary(mod)





### 4. Format NPC data for JAGS run
##subset data to those with NPC known
#npc.data <- bx.full[bx.full$bx.here==1 & !is.na(bx.full$bx.here)
#                    & !is.na(bx.full$npc) & bx.full$time.int>0,]

##number of NPC outcomes
#n_npc <- dim(npc.data)[1]

##NPC outcomes
#NPC <- npc.data$npc

##unique pt id for each bx result
#subj_npc <- npc.data$subj

##covariate matrix
#NCS.offset <- npc.data$ncs
#V.NPC.data <- as.matrix(cbind(rep(1,n_npc),npc.data$mri, npc.data$std.vol))
#d.V.NPC <- dim(V.NPC.data)[2]

##for intervals with second biopsies
#npc.data2 <- bx.full[bx.full$bx.here==1 & !is.na(bx.full$bx.here)
#                     & bx.full$num.bx>1 & !is.na(bx.full$npc.min)
#                     & bx.full$time.int>0,]
#n_npc2 <- dim(npc.data2)[1]
#NPC2 <- npc.data2$npc.min
#subj_npc2 <- npc.data2$subj

#NCS2.offset <- npc.data2$ncs.min
#V.NPC2.data <- as.matrix(cbind(rep(1,n_npc2), npc.data2$mri.min, npc.data2$std.vol))

##combine
#n_npc <- n_npc + n_npc2
#NPC <- c(NPC, NPC2)
#subj_npc <- c(subj_npc, subj_npc2)
#data.check(condition=as.logical(sum(is.na(NPC))==0), message="Missing NPC data. Email Yates; she will check orginal script.")

#NCS.offset<- c(NCS.offset, NCS2.offset)
#npc.ncs.knots <- attr(ns(NCS.offset, 2), "knots")
#npc.ncs.bknots <- attr(ns(NCS.offset, 2), "Boundary.knots")

#V.NPC.data <- rbind(V.NPC.data, V.NPC2.data)
#V.NPC.data <- as.matrix(cbind( ns(NCS.offset,2), V.NPC.data))
#d.V.NPC <- dim(V.NPC.data)[2]
#data.check(condition=as.logical(sum(is.na(V.NPC.data))==0), message="Missing NPC data. Email Yates; she will check orginal script.")


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
n_pgg <- n_pgg + n_pgg2
PGG <- c(PGG, PGG2)
subj_pgg <- c(subj_pgg, subj_pgg2)
data.check(condition=as.logical(sum(is.na(PGG))==0), message="Missing biopsy results data. Email Yates; she will check orginal script.")

V.PGG.data <- rbind(V.PGG.data, V.PGG2.data)
#natural splines for continuous variables
bx.date.knots <- attr(ns(V.PGG.data[,1], 3), "knots")
bx.date.bknots <- attr(ns(V.PGG.data[,1], 3), "Boundary.knots")
pgg.ncs.knots <- attr(ns(V.PGG.data[,2], 2), "knots")
pgg.ncs.bknots <- attr(ns(V.PGG.data[,2], 2), "Boundary.knots")

V.PGG.data <- as.matrix(cbind( ns(V.PGG.data[,1], 3),
                               ns(V.PGG.data[,2], 2),
                               V.PGG.data[,3:d.V.PGG] ))
d.V.PGG <- dim(V.PGG.data)[2]
data.check(condition=as.logical(sum(is.na(V.PGG.data))==0), message="Missing biopsy data. Email Yates; she will check orginal script.")




### 6. Format MPC data for JAGS run
##subset data to biopsies with any cancer found (NPC>0)
#mpc.data <- bx.full[bx.full$bx.here==1 & !is.na(bx.full$bx.here)
#                    & !is.na(bx.full$mpc) & bx.full$mpc>0
#                    & bx.full$time.int>0,]
##number of MPC observations
#n_mpc <- dim(mpc.data)[1]
##MPC findings
#MPC <- mpc.data$mpc
##unique pt identifier
#subj_mpc <- mpc.data$subj

##covariate matrix
#V.MPC.data <- as.matrix(cbind( mpc.data$ncs, mpc.data$mri, mpc.data$std.vol))
#d.V.MPC <- dim(V.MPC.data)[2]

##for intervals with second biopsies with NPC>0
#mpc.data2 <- bx.full[bx.full$bx.here==1 & !is.na(bx.full$bx.here)
#                     & !is.na(bx.full$mpc.min) & bx.full$mpc.min>0
#                     & bx.full$time.int>0,]
#n_mpc2 <- dim(mpc.data2)[1]
#MPC2 <- mpc.data2$mpc.min
#subj_mpc2 <- mpc.data2$subj
#V.MPC2.data <- as.matrix(cbind( mpc.data2$ncs.min, mpc.data2$mri.min, mpc.data2$std.vol))

##combine
#n_mpc <- n_mpc + n_mpc2
#MPC <- c(MPC, MPC2)
#MPC <- round((MPC/10)+0.01)
#subj_mpc <- c(subj_mpc, subj_mpc2)
#data.check(condition=as.logical(sum(is.na(MPC))==0), message="Missing MPC data. Email Yates; she will check orginal script.")

#V.MPC.data <- rbind(V.MPC.data, V.MPC2.data)
#mpc.ncs.knots <- attr(ns(V.MPC.data[,1], 2), "knots")
#mpc.ncs.bknots <- attr(ns(V.MPC.data[,1], 2), "Boundary.knots")

#V.MPC.data <- as.matrix(cbind(rep(1, n_mpc), ns(V.MPC.data[,1],2), V.MPC.data[,2:3] ))
#d.V.MPC <- dim(V.MPC.data)[2]
#data.check(condition=as.logical(sum(is.na(V.MPC.data))==0), message="Missing MPC data. Email Yates; she will check orginal script.")



### 7. Format laterality data for JAGS run
##only include biopsies with NPC>0
#bx.full$lat[bx.full$npc==0]<-NA
#bx.full$lat.min[bx.full$npc.min==0]<-NA
#lat.data <- bx.full[bx.full$bx.here==1 & !is.na(bx.full$bx.here)
#                    & !is.na(bx.full$lat),]
#                    #& bx.full$time.int>0,]
##number observations
#n_lat <- dim(lat.data)[1]
##laterality outcome
#LAT <- lat.data$lat
## unique identifier
#subj_lat <- lat.data$subj

##covariate matrix
#V.LAT.data <- as.matrix(cbind( lat.data$ncs, lat.data$mri, lat.data$std.vol))
#d.V.LAT <- dim(V.LAT.data)[2]

##for intervals with second biopsies
#lat.data2 <- bx.full[bx.full$bx.here==1 & !is.na(bx.full$bx.here)
#                     & !is.na(bx.full$lat.min)
#                     & bx.full$time.int>0,]
#n_lat2 <- dim(lat.data2)[1]
#LAT2 <- lat.data2$lat.min
#subj_lat2 <- lat.data2$subj
#V.LAT2.data <- as.matrix(cbind( lat.data2$ncs.min, lat.data2$mri.min, lat.data2$std.vol))

#combine
#n_lat <- n_lat + n_lat2
#LAT <- c(LAT, LAT2)
#subj_lat <- c(subj_lat, subj_lat2)
#data.check(condition=as.logical(sum(is.na(LAT))==0), message="Missing LAT data. Email Yates; she will check orginal script.")

#V.LAT.data <- rbind(V.LAT.data, V.LAT2.data)
#lat.ncs.knots <- attr(ns(V.LAT.data[,1], 3), "knots")
#lat.ncs.bknots <- attr(ns(V.LAT.data[,1], 3), "Boundary.knots")

#V.LAT.data <- as.matrix(cbind(rep(1, n_lat), ns(V.LAT.data[,1],3), V.LAT.data[,2:3] ))
#d.V.LAT <- dim(V.LAT.data)[2]
#data.check(condition=as.logical(sum(is.na(V.LAT.data))==0), message="Missing LAT data. Email Yates; she will check orginal script.")



### 8. Format TREATMENT received data for JAGS run
#### The decision to get SURGERY or any TREATMENT were considered separately for model inclusion. Neither improved model estimation.
##subset data to not include first bx. we know all patients had at least one post-dc bx
#surg.data <- bx.full[bx.full$time.int>0,]
##binary indicator of decision to get surgery/treatment
#SURG <- surg.data$surg
##number surgery decisions
#n_surg <- dim(surg.data)[1]
## unique identifier
#subj_surg <- surg.data$subj
#data.check(condition=as.logical(sum(is.na(SURG))==0), message="Missing surgery data. Email Yates; she will check orginal script.")

##fill in missing data
#surg.data$prop.npc0.end[is.na(surg.data$prop.npc0.end)]<-0

##covariate matrix
#W.SURG.data <- as.matrix(cbind( rep(1,n_surg),
#                                ns(surg.data$time.int, 4) ,
#                                ns(surg.data$int.dt.num, 3) ,
#                                ns(surg.data$int.age, 3) ,
#                                ns(surg.data$num.prev.bx.end, 2) ,
#                                ns(surg.data$max.prev.mpc.end, 2) ,
#                                surg.data$psa.pres.bx ,
#                                surg.data$std.vol ,
#                                ns(surg.data$prop.npc0.end,2) ,
#                                ns(surg.data$freq.bx.end,2) ,
#                                surg.data$prev.vol.rc.end ,
#                                surg.data$prev.rc ,
#                                surg.data$prev.grade.rc.3 ) )
#d.W.SURG <- dim(W.SURG.data)[2]
#data.check(condition=as.logical(sum(is.na(W.SURG.data))==0), message="Missing surgery covariate data. Email Yates; she will check orginal script.")


#mod<-glm(SURG ~ ns(surg.data$time.int, 4) +
#           ns(surg.data$int.dt.num, 3) +
#           ns(surg.data$int.age, 3) +
#           ns(surg.data$num.prev.bx.end, 2) +
#           ns(surg.data$max.prev.mpc.end, 2) +
##           ns(surg.data$max.prev.npc.end, 2) +
#           surg.data$psa.pres.bx + #surg.data$psa.traj.bx +
#           surg.data$std.vol +
#           surg.data$prev.vol.rc.end +
#           ns(surg.data$prop.npc0.end,2) +
#           ns(surg.data$freq.bx.end,2) +
#           surg.data$prev.rc +
#           surg.data$prev.grade.rc.3 #+
##           #surg.data$prev.grade.rc.4
##          surg.data$prev.lat.end
#           ,
#         family="binomial")
#summary(mod)
