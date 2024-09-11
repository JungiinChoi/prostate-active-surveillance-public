# Computation of AUC using data in generated-files folder


### Define directories, file names
base.location <- workdir
location.of.data <- paste0(base.location, "/data")
location.of.r.scripts <- paste0(base.location, "/code")
location.of.generated.files <- paste0(base.location, "/generated-files")
location.of.generated.folder = paste(location.of.generated.files, "/", mri_role, sep="")

### Load libraries, helper functions
library("pROC")
library("ROCR")
library("splines")

J <- 3

p_eta_hat <- list(length = J)
for(j in 1:J){
  p_eta_hat[[j]] <- t(matrix(apply(read.csv(paste(location.of.generated.folder, "/jags-prediction-p_eta_", 
                                                  j,"-", mri_role,".csv",sep=""))[,-1], 2, mean), ncol = 4))
}

# AUC dataframe for PGG > 1, PGG > 2, PGG > 3
# Plot ROC Curve

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

## AUC results
AUC_list

