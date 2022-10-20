get_stats <- function(x){a <- quantile(x, c(0.025, 0.975)); b<-mean(x); return(c(a[1], b, a[2]))}
expit <- function(x){exp(x)/(1+exp(x))}
library(readr)
# pirads_int1 <-read_csv("generated-files/generated-files-new/seq-mri-outcome/jags-prediction-pirads_int1-2022.csv")
# pirads_int2 <-read_csv("generated-files/generated-files-new/seq-mri-outcome/jags-prediction-pirads_int2-2022.csv")
# pirads_slope1 <-read_csv("generated-files/generated-files-new/seq-mri-outcome/jags-prediction-pirads_slope1-2022.csv")
# pirads_slope2 <-read_csv("generated-files/generated-files-new/seq-mri-outcome/jags-prediction-pirads_slope2-2022.csv")
# 
# 
# apply(pirads_int1, 2, get_stats)
# apply(pirads_int2, 2, get_stats)
# apply(pirads_slope1, 2, get_stats)
# apply(pirads_slope2, 2, get_stats)

pgg_slope1 <-  read_csv("generated-files/generated-files-new/seq-mri-moderator/jags-prediction-pgg_slope1-2022.csv")
pgg_slope2 <-  read_csv("generated-files/generated-files-new/seq-mri-moderator/jags-prediction-pgg_slope2-2022.csv")
pgg_slope3 <-  read_csv("generated-files/generated-files-new/seq-mri-moderator/jags-prediction-pgg_slope3-2022.csv")

apply(pgg_slope1, 2, get_stats)
apply(pgg_slope2, 2, get_stats)
apply(pgg_slope3, 2, get_stats)

pgg_slope1 <-  read_csv("generated-files/generated-files-new/seq-mri-moderator-form3/jags-prediction-pgg_slope1-2022.csv")
pgg_slope2 <-  read_csv("generated-files/generated-files-new/seq-mri-moderator-form3/jags-prediction-pgg_slope2-2022.csv")
pgg_slope3 <-  read_csv("generated-files/generated-files-new/seq-mri-moderator-form3/jags-prediction-pgg_slope3-2022.csv")

apply(pgg_slope1, 2, get_stats)
apply(pgg_slope2, 2, get_stats)
apply(pgg_slope3, 2, get_stats)


pgg_slope1 <-  read_csv("generated-files/generated-files-new/seq-mri-both/jags-prediction-pgg_slope1-2022.csv")
pgg_slope2 <-  read_csv("generated-files/generated-files-new/seq-mri-both/jags-prediction-pgg_slope2-2022.csv")
pgg_slope3 <-  read_csv("generated-files/generated-files-new/seq-mri-both/jags-prediction-pgg_slope3-2022.csv")

apply(pgg_slope1, 2, get_stats)
apply(pgg_slope2, 2, get_stats)
apply(pgg_slope3, 2, get_stats)

pgg_slope1 <-  read_csv("generated-files/generated-files-new/seq-mri-both-form3/jags-prediction-pgg_slope1-2022.csv")
pgg_slope2 <-  read_csv("generated-files/generated-files-new/seq-mri-both-form3/jags-prediction-pgg_slope2-2022.csv")
pgg_slope3 <-  read_csv("generated-files/generated-files-new/seq-mri-both-form3/jags-prediction-pgg_slope3-2022.csv")

apply(pgg_slope1, 2, get_stats)
apply(pgg_slope2, 2, get_stats)
apply(pgg_slope3, 2, get_stats)

files_pgs1 <- list.files(path="generated-files-cv/new/seq-mri-moderator", pattern="pgg_slope3-", full.names = T)
files_pgs1_list <- lapply(files_pgs1, read_csv)
pgs1_summ_list <- lapply(files_pgs1_list, apply, 2, get_stats)
pgs1_summ <- NULL
for(i in 1:length(pgs1_summ_list)){
  tmp <- pgs1_summ_list[[i]][,"V11"]
  pgs1_summ <- rbind(pgs1_summ, tmp)
}
apply(pgs1_summ,2, mean)