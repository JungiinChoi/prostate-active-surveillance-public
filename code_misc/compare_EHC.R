## check significant parameter in C bgg model 
library(readr)
c_bgg_int1 <- read_csv("generated-files-seq/jags-prediction-pgg_int1-2022.csv")
c_bgg_int2 <- read_csv("generated-files-seq/jags-prediction-pgg_int2-2022.csv")
c_bgg_int3 <- read_csv("generated-files-seq/jags-prediction-pgg_int3-2022.csv")

c_bgg_slope1 <- read_csv("generated-files-seq/jags-prediction-pgg_slope1-2022.csv")
c_bgg_slope2 <- read_csv("generated-files-seq/jags-prediction-pgg_slope2-2022.csv")
c_bgg_slope3 <- read_csv("generated-files-seq/jags-prediction-pgg_slope3-2022.csv")

c_bgg_int1 <- c_bgg_int1[,-1]
c_bgg_int2 <- c_bgg_int2[,-1]
c_bgg_int3 <- c_bgg_int3[,-1]
c_bgg_slope1 <- c_bgg_slope1[,-1]
c_bgg_slope2 <- c_bgg_slope2[,-1]
c_bgg_slope3 <- c_bgg_slope3[,-1]

get_stats <- function(x){a <- quantile(x, c(0.025, 0.975)); b<-mean(x); return(c(a[1], b, a[2]))}
apply(c_bgg_int1, 2, get_stats)
apply(c_bgg_int2, 2, get_stats)
apply(c_bgg_int3, 2, get_stats)

apply(c_bgg_slope1, 2, get_stats)
apply(c_bgg_slope2, 2, get_stats)
apply(c_bgg_slope3, 2, get_stats)

## spline plot for C showing how biopsy changes w/ time
source("R-script/data-prep-for-jags-reform-seq-mri-abserror.R")
bgg_int1 <- apply(c_bgg_int1, 2, mean)
bgg_int2 <- apply(c_bgg_int2, 2, mean)
bgg_int3 <- apply(c_bgg_int3, 2, mean)
bgg_slope1 <- apply(c_bgg_slope1, 2, mean)
bgg_slope2 <- apply(c_bgg_slope2, 2, mean)
bgg_slope3 <- apply(c_bgg_slope3, 2, mean)

cancer_state <- read_csv("generated-files-seq/jags-prediction-cancer_state-2022.csv")
expit <- function(x) exp(x)(1+exp(x))

## E and H
e_pirads_int1 <- read_csv("generated-files-seq-pgg-removed/jags-prediction-pirads_int1-2022.csv")
e_pirads_int2 <- read_csv("generated-files-seq-pgg-removed/jags-prediction-pirads_int2-2022.csv")
e_pirads_slope1 <- read_csv("generated-files-seq-pgg-removed/jags-prediction-pirads_slope1-2022.csv")
e_pirads_slope2 <- read_csv("generated-files-seq-pgg-removed/jags-prediction-pirads_slope2-2022.csv")
apply(e_pirads_int1, 2, get_stats)
apply(e_pirads_int2, 2, get_stats)
apply(e_pirads_slope1, 2, get_stats)
apply(e_pirads_slope2, 2, get_stats)


h_pirads_int1 <- read_csv("generated-files-seq-mri-abserror/jags-prediction-pirads_int1-2022.csv")
h_pirads_int2 <- read_csv("generated-files-seq-mri-abserror/jags-prediction-pirads_int2-2022.csv")
h_pirads_slope1 <- read_csv("generated-files-seq-mri-abserror/jags-prediction-pirads_slope1-2022.csv")
h_pirads_slope2 <- read_csv("generated-files-seq-mri-abserror/jags-prediction-pirads_slope2-2022.csv")
apply(h_pirads_int1, 2, get_stats)
apply(h_pirads_int2, 2, get_stats)
apply(h_pirads_slope1, 2, get_stats)
apply(h_pirads_slope2, 2, get_stats)

#H model abserror term
h_abserror_int1 <- read_csv("generated-files-seq-mri-abserror/jags-prediction-abserror_int1-2022.csv")
h_abserror_int2 <- read_csv("generated-files-seq-mri-abserror/jags-prediction-abserror_int2-2022.csv")
h_abserror_int3 <- read_csv("generated-files-seq-mri-abserror/jags-prediction-abserror_int3-2022.csv")
h_abserror_slope <- read_csv("generated-files-seq-mri-abserror/jags-prediction-abserror_slope-2022.csv")

apply(h_abserror_int1, 2, get_stats)
apply(h_abserror_int2, 2, get_stats)
apply(h_abserror_int3, 2, get_stats)
apply(h_abserror_slope, 2, get_stats)
