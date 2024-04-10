# Computation of AUC using data in generated-files folder


### Clear workspace
rm(list=ls())


### Load libraries, helper functions
library("pROC")
library("ROCR")
library("splines")

p_eta_hat <- readRDS(paste0(location.of.generated.folder, "/p_eta_hat.rds"))

# AUC dataframe for PGG > 1, PGG > 2, PGG > 3
# Plot ROC Curve

J <- 3
dx_list <- list(length = J)
for (i in 1:J){
  dx_list[[i]] <- read.csv(paste0(location.of.data,"/dx_",i,".csv"))}

par(mfrow = c(J,3))
AUC_list <- data.frame(dx_1 = rep(0,3), dx_2 = rep(0,3), dx_3 = rep(0,3))
for (i in 1:J){
  ind <- is.na(dx_list[[i]]$rp)
  obs <- dx_list[[i]]$rp[!ind]
  obs_larger <- cbind(as.numeric(obs>1), as.numeric(obs>2), as.numeric(obs>3))
  exp <- p_eta_hat[[i]][, !ind]
  exp_larger <- t(apply(exp, 2, function(x){c(sum(x[2:4]), sum(x[3:4]), x[4])}))
  for (k in 1:3){
    roc_tmp <- roc(obs_larger[,k], exp_larger[,k])
    AUC_list[k,i] <- auc(roc_tmp)
    plot(roc_tmp, main = paste0("Cohort ", i, ": PGG > ", k))
  }
}

