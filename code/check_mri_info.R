setwd("~/Downloads/prostate-active-surveillance-vDaan")
options(warn=1)
data.check <- function(condition, message){
  if(condition==FALSE){print(paste(message, "Program terminated.", sep=" "))}
  stopifnot(condition)}
load("~/Downloads/prostate-active-surveillance-vDaan/generated-files-cv/IOP-data-shaping-work-space.RData")
source("R-script/data-prep-for-jags-reform-addmri.R")

length(unique(mri_data$clinical_ptnum))
length(unique(pt.data$clinical_PTnum))
sum(unique(mri_data$clinical_ptnum) %in% unique(pt.data$clinical_PTnum))
id_etaknown <- unique(pt.data$clinical_PTnum[!is.na(pt.data$true.pgg)])
sum(unique(mri_data$clinical_ptnum) %in% id_etaknown)

mri_data$pirads_max
tmp_dt <- pt.data %>% 
  dplyr::select(clinical_PTnum, true.pgg, bx.pgg) %>% 
  left_join(mri_data, by = c("clinical_PTnum" = "clinical_ptnum"))
tmp_dt <- tmp_dt %>% filter(!is.na(true.pgg))
tmp_dt$diff_eta_pgg <- abs(tmp_dt$true.pgg - tmp_dt$bx.pgg)
tmp_dt$have_mri <- ifelse(is.na(tmp_dt$pirads_max), 'no', 'yes')
t <- with(tmp_dt, table(diff_eta_pgg, have_mri)) %>%
  prop.table(margin = 2)
round(t,2)
table(tmp_dt$diff_eta_pgg, tmp_dt$have_mri)

table(tmp_dt$diff_eta_pgg, tmp_dt$pirads_max)
round(with(tmp_dt, table(diff_eta_pgg, pirads_max)) %>%
  prop.table(margin = 2),2)

pirads_data <- tmp_dt$pirads_max[!is.na(tmp_dt$pirads_max)]
sum(pirads_data %in% c(1,2))/length(pirads_data)
sum(pirads_data %in% c(3))/length(pirads_data)
sum(pirads_data %in% c(4,5))/length(pirads_data)
logit <- function(p) log(p/(1-p))

logit(sum(pirads_data %in% c(1,2))/length(pirads_data))
logit(sum(pirads_data %in% c(3))/sum(pirads_data %in% c(3,4,5)))

##
length(tmp_dt$diff_eta_pgg[!is.na(tmp_dt$diff_eta_pgg)])
table(tmp_dt$diff_eta_pgg)
abserror <- tmp_dt$diff_eta_pgg

logit(sum(abserror == 0)/length(abserror))
logit(sum(abserror == 1)/sum(abserror %in% c(1,2,3)))
logit(sum(abserror == 2)/sum(abserror %in% c(2,3)))

##
tmp_layer2 <- tmp_dt %>% 
  filter(diff_eta_pgg == 1 & true.pgg == 2) %>% 
  mutate(eta_ge_pgg = as.numeric(true.pgg > bx.pgg))
tmp_layer3 <- tmp_dt %>% 
  filter(diff_eta_pgg == 1 & true.pgg == 3) %>% 
  mutate(eta_ge_pgg = as.numeric(true.pgg > bx.pgg))

layer_eta2_data <- tmp_layer2$eta_ge_pgg
layer_eta3_data <- tmp_layer3$eta_ge_pgg

logit(sum(layer_eta2_data == 1)/length(layer_eta2_data))
logit(sum(layer_eta3_data == 1)/length(layer_eta3_data))


