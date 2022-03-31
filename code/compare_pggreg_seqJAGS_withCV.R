## compare sequential model and proportional odds using patients under surgery (eta known)
## for both varying prior and fixed intercept variance
## added cross validation
rm(list=ls())
source('R-scripts/load_libs.R')
base.location <- "/Users/zitongwang/Downloads/prostate-active-surveillance-vDaan/" #"/users/rcoley/jhas-epic/psa/psa-new/" #"/Users/ryc/Documents/inhealth/prediction-model/automated/for-TIC/"
setwd(base.location)
location.of.data <- paste0(base.location, "data")
location.of.r.scripts <- paste0(base.location, "R-scripts")

# name.of.pt.data <- "Demographic_6.15.csv" #"demographics with physician info.2015.csv"
# name.of.bx.data <- "Biopsy_6.15.csv" #"Biopsy data_2015.csv"
# name.of.psa.data <- "PSA_6.15.csv" #"PSA.2015.csv"
# name.of.tx.data <- "Treatment_6.15.csv" #"treatment_2015.csv"
# date.pull<-as.numeric(Sys.Date())
# source(paste(location.of.r.scripts,"data-load-check-and-shaping.R",sep="/"))
# save.image(file='dataloading.RData')
load('dataloading.RData')
location.of.generated.files <- paste0(base.location, "generated-files-pgg-model/sequential_model/cv")
#model.file = paste("code/compare_pggreg_sigma_seqJAGS.txt")
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
# pgg.data$eta_m2 <- ifelse(pgg.data$true.pgg-2 >=0 , 1, 0)
# pgg.data$eta_m3 <- ifelse(pgg.data$true.pgg-3 >=0 , 1, 0)
# pgg.data$eta_m4 <- ifelse(pgg.data$true.pgg-4 >=0 , 1, 0)

bx.date.knots <- attr(ns(pgg.data$bx.dt.num, 3), "knots")
bx.date.bknots <- attr(ns(pgg.data$bx.dt.num, 3), "Boundary.knots")
pgg.ncs.knots <- attr(ns(pgg.data$ncs, 2), "knots")
pgg.ncs.bknots <- attr(ns(pgg.data$ncs, 2), "Boundary.knots")

std.ns.bxdt <- scale(ns(pgg.data$bx.dt.num, Boundary.knots = bx.date.bknots, knots = bx.date.knots))
std.ns.ncs <- scale(ns(pgg.data$ncs , Boundary.knots = pgg.ncs.bknots, knots = pgg.ncs.knots))
modmat_pgg <- as.matrix(cbind(std.ns.bxdt,
                              std.ns.ncs,
                              pgg.data[, c("mri", "std.vol")] ))
npred_pgg <- dim(modmat_pgg)[2]
colnames(modmat_pgg) <- c("ns(bx.dt.num,1)","ns(bx.dt.num,2)","ns(bx.dt.num,3)",
                          "ns(std.ncs,1)", "ns(std.ncs,2)", 
                          "mri", "std.vol")
pgg_data <- pgg.data$pgg
npat_pgg <- n_pgg1 + n_pgg2
#The number of latent classes/ values of true cancer state
nlevel_cancer <- 4

## cross validate
idlist1 <- unique(pgg.data$clinical_PTnum[pgg.data$pgg==1])
idlist2 <- unique(pgg.data$clinical_PTnum[pgg.data$pgg==2])
idlist3 <- unique(pgg.data$clinical_PTnum[pgg.data$pgg==3])
idlist4 <- unique(pgg.data$clinical_PTnum[pgg.data$pgg==4])
set.seed(2022)
test_inx1 <- createFolds(idlist1, k=5, list = T)
test_inx2 <- createFolds(idlist2, k=5, list = T)
test_inx3 <- createFolds(idlist3, k=5, list = T)
test_inx4 <- createFolds(idlist4, k=5, list = T)

test_ids1 <- test_ids2 <- test_ids3 <- test_ids4 <- list()
get_testid <- function(testinx, idlist){
  test_ids <- list()
  for(i in 1:length(testinx)){
    test.inx <- testinx[[i]]
    test_ids[[i]] <- idlist[test.inx]
  }
  return(test_ids)
}
test_ids1 <- get_testid(test_inx1, idlist1)
test_ids2 <- get_testid(test_inx2, idlist2)
test_ids3 <- get_testid(test_inx3, idlist3)
test_ids4 <- get_testid(test_inx4, idlist4)

### run model
coef_mean_opt <-"betaj" #betaj(k-1), betaj
coef_var_opt <- "sigma"  #sigma, sigmaj
model.file = paste0("code/jags_model/pggreg_seqJAGS_",
                    coef_mean_opt, "_", coef_var_opt,".txt")
location.of.generated.files <- 
  paste0(base.location, 
         "generated-files-pgg-model/sequential_model/cv/", 
         coef_mean_opt, "_", coef_var_opt)

for(k in 1:5){
  message(k)
  test_id_stage1 = c(test_ids1[[k]],test_ids2[[k]],test_ids3[[k]],test_ids4[[k]])
  train_pos1 = which(!(pgg.data$clinical_PTnum[pgg.data$pgg>=1] %in% test_id_stage1))
  pgg_data1 = pgg_data[train_pos1]; pgg_data1 = ifelse(pgg_data1 == 1, 1, 0)
  modmat_pgg1 = modmat_pgg[train_pos1,]
  npat_pgg1 = length(pgg_data1)
  
  test_id_stage2 = c(test_ids2[[k]],test_ids3[[k]],test_ids4[[k]])
  train_pos2 = which(!(pgg.data$clinical_PTnum[pgg.data$pgg>=2] %in% test_id_stage2))
  inx_lev2 <- which(pgg_data > 1)
  pgg_data2 = pgg_data[inx_lev2]; pgg_data2 = pgg_data2[train_pos2]
  pgg_data2 = ifelse(pgg_data2 == 2, 1, 0)
  modmat_pgg2 = modmat_pgg[inx_lev2,]; modmat_pgg2 = modmat_pgg2[train_pos2,]
  npat_pgg2 = length(pgg_data2)
  
  test_id_stage3 = c(test_ids3[[k]],test_ids4[[k]])
  train_pos3 = which(!(pgg.data$clinical_PTnum[pgg.data$pgg>=3] %in% test_id_stage3))
  inx_lev3 <- which(pgg_data > 2)
  pgg_data3 = pgg_data[inx_lev3]; pgg_data3 = pgg_data3[train_pos3]
  pgg_data3 = ifelse(pgg_data3 == 3, 1, 0)
  modmat_pgg3 = modmat_pgg[inx_lev3,]; modmat_pgg3 = modmat_pgg3[train_pos3,]
  npat_pgg3 = length(pgg_data3)

  ### 1. Set up jags arguments
  jags_data<-list(
    pgg_data1=pgg_data1, pgg_data2=pgg_data2,pgg_data3=pgg_data3,  
    modmat_pgg1=modmat_pgg1, modmat_pgg2=modmat_pgg2,modmat_pgg3=modmat_pgg3,  
    npat_pgg1=npat_pgg1,npat_pgg2=npat_pgg2,npat_pgg3=npat_pgg3, 
    npred_pgg = npred_pgg,
    alpha = 0.001, beta = 0.001
  ) 
  ### 2. Initialize model parameters
  inits <- function() {
    pgg_int1 <- pgg_int2 <-pgg_int3 <-rnorm(1,0,1)
    pgg_slope1 <- rnorm(npred_pgg,mean=0,sd=0.25)
    pgg_slope2 <- rnorm(npred_pgg,mean=0,sd=0.25)
    pgg_slope3 <- rnorm(npred_pgg,mean=0,sd=0.25)
    if(coef_var_opt == "sigma"){
      inv_sigma_sq = rgamma(1, 1, 1)
    }else if(coef_var_opt == "sigmaj"){
      inv_sigma_sq = rgamma(npred_pgg, 1, 1)
    }
    if(coef_mean_opt == "betaj"){
      coef_mean = rnorm(npred_pgg,mean=0,sd=0.25)
      list(
        pgg_int1=pgg_int1, pgg_int2=pgg_int2, pgg_int3=pgg_int3, 
        pgg_slope1=pgg_slope1,pgg_slope2=pgg_slope2,pgg_slope3=pgg_slope3,
        inv_sigma_sq = inv_sigma_sq, coef_mean = coef_mean
      )
    }else{
      list(
        pgg_int1=pgg_int1, pgg_int2=pgg_int2, pgg_int3=pgg_int3, 
        pgg_slope1=pgg_slope1,pgg_slope2=pgg_slope2,pgg_slope3=pgg_slope3,
        inv_sigma_sq = inv_sigma_sq
      )
    }
    
     }
  
    ### 3. Define parameters to be tracked
    params <- c(
      "pgg_int1", "pgg_int2", "pgg_int3",
      "pgg_slope1", "pgg_slope2", "pgg_slope3", "inv_sigma_sq",
      "sigma_sq"
    )

  
  ### 4. Define other jags settings  
  n.iter <- 10000; n.burnin <- 2500; n.thin <- 10; n.chains <- 1
  
  ### 5. Run
  seed=2022
  set.seed(seed)
  outj<-jags(jags_data, inits=inits, parameters.to.save=params,
             model.file=model.file,
             n.thin=n.thin, n.chains=n.chains, n.burnin=n.burnin, n.iter=n.iter)
  out<-outj$BUGSoutput
  
  for(j in 1:length(out$sims.list)){
    write.csv(out$sims.list[[j]],
              paste(location.of.generated.files, 
                    "/jags_", names(out$sims.list)[j],"_",seed,"cv",k,".csv",sep=""))
    
  }
}

## average output
get_stats <- function(x){a <- quantile(x, c(0.025, 0.975)); b<-mean(x); return(c(a[1], b, a[2]))}
expit <- function(x){exp(x)/(1+exp(x))}
nsim <- ((n.iter - n.burnin)/n.thin)
exp_cv <- obs_cv <- list()
coef_mean_opt <- "betaj" #betaj, betaj(k-1)
coef_var_opt <- "sigma"  #sigma, sigmaj
location.of.generated.files <-  paste0(base.location, 
                                       "generated-files-pgg-model/sequential_model/cv/", 
                                       coef_mean_opt, "_", coef_var_opt)
for(k in 1:5){
  jags_int1 <-
    read_csv(paste(location.of.generated.files, "/jags_pgg_int1_",seed,"cv",k,".csv",sep=""))
  jags_int2 <-
    read_csv(paste(location.of.generated.files, "/jags_pgg_int2_",seed,"cv",k,".csv",sep=""))
  jags_int3 <-
    read_csv(paste(location.of.generated.files, "/jags_pgg_int3_",seed,"cv",k,".csv",sep=""))
  jags_slope1 <-
    read_csv(paste(location.of.generated.files, "/jags_pgg_slope1_",seed,"cv",k,".csv",sep=""))
  jags_slope2 <-
    read_csv(paste(location.of.generated.files, "/jags_pgg_slope2_",seed,"cv",k,".csv",sep=""))
  jags_slope3 <-
    read_csv(paste(location.of.generated.files, "/jags_pgg_slope3_",seed,"cv",k,".csv",sep=""))
  jags_int1 <- as.matrix(jags_int1[,-1])
  jags_int2 <- as.matrix(jags_int2[,-1])
  jags_int3 <- as.matrix(jags_int3[,-1])
  jags_slope1 <- as.matrix(jags_slope1[,-1])
  jags_slope2 <- as.matrix(jags_slope2[,-1])
  jags_slope3 <- as.matrix(jags_slope3[,-1])
  ## predict on test sets
  test_id_all = c(test_ids1[[k]],test_ids2[[k]],test_ids3[[k]],test_ids4[[k]])
  test_pos = which(pgg.data$clinical_PTnum %in% test_id_all)
  modmat_pgg_test = modmat_pgg[test_pos,]
  pgg_data_test = pgg.data[test_pos,]

  p_rc1 <-p_rc2 <-p_rc3 <- p_rc4 <-matrix(0, ncol = nsim, nrow = dim(modmat_pgg_test)[1])
  for(i in 1:nsim){
    int1 <- matrix(jags_int1[i,])
    int2 <-  matrix(jags_int2[i,])
    int3 <-  matrix(jags_int3[i,])
    slope1 <- matrix(jags_slope1[i,])
    slope2 <- matrix(jags_slope2[i,])
    slope3 <- matrix(jags_slope3[i,])

    linpred1 <- modmat_pgg_test %*% slope1
    p_y1_given_x0 <- expit(matrix(int1 + c(linpred1)))
    linpred2 <- modmat_pgg_test %*% slope2
    p_y2_given_x0ge1 <- expit(matrix(int2 + c(linpred2)))
    linpred3 <- modmat_pgg_test %*% slope3
    p_y3_given_x0ge2 <- expit(matrix(int3 + c(linpred3)))

    p_rc1[,i] <- p_y1_given_x0
    p_rc2[,i] <- (1-p_y1_given_x0) * p_y2_given_x0ge1
    p_rc3[,i] <- (1-p_y1_given_x0)*(1-p_y2_given_x0ge1) * p_y3_given_x0ge2
    p_rc4[,i] <- 1- p_rc1[,i] -  p_rc2[,i] -  p_rc3[,i]
  }
  exp1 <- sum(apply(p_rc1, 1, mean))
  exp2 <- sum(apply(p_rc2, 1, mean))
  exp3 <- sum(apply(p_rc3, 1, mean))
  exp4 <- sum(apply(p_rc4, 1, mean))
  obs1 <- sum(pgg_data_test$pgg == 1);obs2 <- sum(pgg_data_test$pgg == 2)
  obs3 <- sum(pgg_data_test$pgg == 3);obs4 <- sum(pgg_data_test$pgg == 4)
  exp_cv[[k]] <- list(exp1=exp1, exp2=exp2, exp3=exp3, exp4=exp4)
  obs_cv[[k]] <- list(obs1=obs1, obs2=obs2, obs3=obs3, obs4=obs4)
}

location.of.generated.files
exp_ave <- apply(apply(do.call(rbind, exp_cv), 2,as.numeric),2, mean)
obs_ave <- apply(apply(do.call(rbind, obs_cv), 2,as.numeric),2, mean)
location.of.generated.files
obs_ave; exp_ave
sum((obs_ave-exp_ave)^2/exp_ave)

# compare variance across models -------------
set.seed(2022); library(invgamma);prior_distr <- rinvgamma(750,1,1)
#get posterior distribution
variable_extract <- "sigma_sq"
coef_seq <- c("betaj(k-1)_sigmaj", "betaj(k-1)_sigma",
              "betaj_sigmaj", "betaj_sigma")
jags_post_allmodel <- list()
for(i in 1:length(coef_seq)){
  coef_type <- coef_seq[i]
    location.of.generated.files <- 
      paste0(base.location, 
             "generated-files-pgg-model/sequential_model/cv/", 
             coef_type)
    input_dir <- location.of.generated.files
    input_files <- list.files(path = input_dir, pattern = variable_extract, full.names = TRUE)
    jags_cvlist <- list()
    for(inx in 1:length(input_files)){
      tmp <- read_csv(input_files[inx])[,-1]
      jags_cvlist[[inx]] <- apply(tmp, 2, as.numeric)
    }
    jags_cv_mean <- Reduce(`+`, jags_cvlist)/length(jags_cvlist)
    jags_post_allmodel[[i]] <- jags_cv_mean
    names(jags_post_allmodel)[i] = coef_type
}

par(mfrow = c(3,3), oma = c(0,0,0,0), mar = c(2,2,2,2))
for(icol in 1:7){
  ymax = range(c(density(sqrt(prior_distr))$y, 
                 density(sqrt(jags_post_allmodel[[1]][, icol]))$y,
                 density(sqrt(jags_post_allmodel[[3]][, icol]))$y))
  xlim = range(c(density(sqrt(jags_post_allmodel[[1]][, icol]))$x,
                 density(sqrt(jags_post_allmodel[[3]][, icol]))$x))
  plot(density(sqrt(prior_distr)), col = "grey", type = "l",
       ylim = c(0,ymax[2]),
       xlim = xlim,
       main = paste0(variable_extract, " coef_", icol), xlab = "")
  lines(density(sqrt(jags_post_allmodel[[1]][, icol])), col = "lightblue")
  lines(density(sqrt(jags_post_allmodel[[3]][, icol])), col = "lightpink")
  legend(x = "topright",          # Position
         legend = c("prior", "betaj(k-1)_sigmaj", "betaj_sigmaj"),  # Legend texts
         lty = c(1, 1, 1),           # Line types
         col = c("grey", "lightblue", "lightpink"),           # Line colors
         lwd = 1,
         cex= 0.8) 
}
!(grepl("sigmaj", names(jags_post_allmodel)))
ymax = range(c(density(sqrt(prior_distr))$y, 
               density(sqrt(jags_post_allmodel[[2]]))$y,
               density(sqrt(jags_post_allmodel[[4]]))$y))
xlim = range(c(density(sqrt(jags_post_allmodel[[2]]))$x,
               density(sqrt(jags_post_allmodel[[4]]))$x))

plot(density(sqrt(prior_distr)), col = 'grey', type = 'l', ylim = c(0, ymax[2]),
      main = paste(variable_extract), xlab = " ")
lines(density(sqrt(jags_post_allmodel[[2]])), col = 'lightblue')
lines(density(sqrt(jags_post_allmodel[[4]])), col = 'lightpink')
legend(x = "topright",          # Position
         legend = c("prior", "betaj(k-1)_sigma", "betaj_sigma"),  # Legend texts
         lty = c(1, 1, 1),           # Line types
         col = c("grey", "lightblue", "lightpink"),           # Line colors
         lwd = 1,
         cex= 0.8) 

# compare 3rd stage intercept across models -------------
set.seed(2022); prior_distr <- rnorm(750,0,1)
#get posterior distribution
variable_extract <- "pgg_int3"
coef_seq <- c("betaj(k-1)_sigmaj", "betaj(k-1)_sigma",
              "betaj_sigmaj", "betaj_sigma")
jags_post_allmodel <- list()
for(i in 1:length(coef_seq)){
  coef_type <- coef_seq[i]
  location.of.generated.files <- 
    paste0(base.location, 
           "generated-files-pgg-model/sequential_model/cv/", 
           coef_type)
  input_dir <- location.of.generated.files
  input_files <- list.files(path = input_dir, pattern = variable_extract, full.names = TRUE)
  jags_cvlist <- list()
  for(inx in 1:length(input_files)){
    tmp <- read_csv(input_files[inx])[,-1]
    jags_cvlist[[inx]] <- apply(tmp, 2, as.numeric)
  }
  jags_cv_mean <- Reduce(`+`, jags_cvlist)/length(jags_cvlist)
  jags_post_allmodel[[i]] <- jags_cv_mean
  names(jags_post_allmodel)[i] = coef_type
}

yrange = range(c(density(prior_distr)$y, 
               density(jags_post_allmodel[[1]])$y,
               density(jags_post_allmodel[[2]])$y,
               density(jags_post_allmodel[[3]])$y,
               density(jags_post_allmodel[[4]])$y))
xrange = range(c(density(prior_distr)$x, 
               density(jags_post_allmodel[[1]])$x,
               density(jags_post_allmodel[[2]])$x,
               density(jags_post_allmodel[[3]])$x,
               density(jags_post_allmodel[[4]])$x))
plot(density(prior_distr), col = 'grey', type = 'l', ylim = yrange,
     xlim = xrange,
     main = paste(variable_extract), xlab = " ", lty = 1)
lines(density(jags_post_allmodel[[1]]), col = 'lightblue', lty = 1)
lines(density(jags_post_allmodel[[2]]), col = 'lightblue', lty = 2)
lines(density(jags_post_allmodel[[3]]), col = 'lightpink', lty = 1)
lines(density(jags_post_allmodel[[4]]), col = 'lightpink', lty = 2)
legend(x = "topleft",          # Position
       legend = c("prior",names(jags_post_allmodel)),  # Legend texts
       lty = c(1,1, 2, 1,2),           # Line types
       col = c("grey", "lightblue", "lightblue","lightpink", "lightpink"),           # Line colors
       lwd = 1,
       cex= 0.8) 

# compare 3rd stage coefficient across models -------------
set.seed(2022); prior_distr <- rnorm(750,0,1)
#get posterior distribution
variable_extract <- "pgg_slope3"
coef_seq <- c("betaj(k-1)_sigmaj", "betaj(k-1)_sigma",
              "betaj_sigmaj", "betaj_sigma")
jags_post_allmodel <- list()
for(i in 1:length(coef_seq)){
  coef_type <- coef_seq[i]
  location.of.generated.files <- 
    paste0(base.location, 
           "generated-files-pgg-model/sequential_model/cv/", 
           coef_type)
  input_dir <- location.of.generated.files
  input_files <- list.files(path = input_dir, pattern = variable_extract, full.names = TRUE)
  jags_cvlist <- list()
  for(inx in 1:length(input_files)){
    tmp <- read_csv(input_files[inx])[,-1]
    jags_cvlist[[inx]] <- apply(tmp, 2, as.numeric)
  }
  jags_cv_mean <- Reduce(`+`, jags_cvlist)/length(jags_cvlist)
  jags_post_allmodel[[i]] <- jags_cv_mean
  names(jags_post_allmodel)[i] = coef_type
}


par(mfrow = c(3,3), oma = c(0,0,0,0), mar = c(2,2,2,2))
for(icol in 1:7){
  yrange = range(c(density(prior_distr)$y, 
                   density(jags_post_allmodel[[1]][, icol])$y,
                   density(jags_post_allmodel[[2]][, icol])$y,
                   density(jags_post_allmodel[[3]][, icol])$y,
                   density(jags_post_allmodel[[4]][, icol])$y))
  xrange = range(c(density(prior_distr)$x, 
                   density(jags_post_allmodel[[1]][, icol])$x,
                   density(jags_post_allmodel[[2]][, icol])$x,
                   density(jags_post_allmodel[[3]][, icol])$x,
                   density(jags_post_allmodel[[4]][, icol])$x))
  plot(density(prior_distr), col = 'grey', type = 'l', ylim = yrange,
       xlim = xrange,
       main = paste(variable_extract), xlab = " ", lty = 1)
  lines(density(jags_post_allmodel[[1]][, icol]), col = 'lightblue', lty = 1)
  lines(density(jags_post_allmodel[[2]][, icol]), col = 'lightblue', lty = 2)
  lines(density(jags_post_allmodel[[3]][, icol]), col = 'lightpink', lty = 1)
  lines(density(jags_post_allmodel[[4]][, icol]), col = 'lightpink', lty = 2)
  legend(x = "topright",          # Position
         legend = c("prior",names(jags_post_allmodel)),  # Legend texts
         lty = c(1,1, 2, 1,2),           # Line types
         col = c("grey", "lightblue", "lightblue","lightpink", "lightpink"),           # Line colors
         lwd = 1,
         cex= 0.8)  
}



## betaj_sigmaj