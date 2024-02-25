# Synthesize Data
library(readxl)
location.of.ucsf.data <- paste(location.of.data, "UCSF_ASRP_COHORT_DATA_01MAY2023_MERGED.xlsx", sep="/")
location.of.ucsf.psa.data <- paste(location.of.data, "UCSF_ASRP_COHORT_DATA_01MAY2023_N_360.xlsx", sep="/")
dx.ucsf <- read_excel(location.of.ucsf.data, sheet = 1)
mri.ucsf <- read_excel(location.of.ucsf.data, sheet = 3)
bxmri.ucsf <- read_excel(location.of.ucsf.data, sheet = 4)
psa.ucsf <- read_excel(location.of.ucsf.psa.data, sheet = 4)
bx.ucsf <- read_excel(location.of.ucsf.data, sheet = 2)

# Generate dx_1, bx_1, bxmri_1, psa_1: mask ID of UCSF data
set.seed(011824)
## dx_1: Remove data after RP

dx_1 <- dx.ucsf
dx_1$ID <- sample(1:nrow(dx_1), nrow(dx_1))
ID_mapping <- dx_1 %>% select(PatientID, ID, dxrpdays)
dx_1 <- dx_1[,c(8, (2:7))]
dx_1 <- dx_1[order(dx_1$ID),]
dx_1$vol <- rep(0, nrow(dx_1))
for (i in 1:nrow(dx_1)){
  dx_1$vol[i] <- rnorm(1, 50 + dx_1$rp[i], 10)
}
dx_1[350:360, 2:4] <- NA

## bx_1

bx_1 <- bx.ucsf
bx_1 <- bx_1 %>% 
  left_join(ID_mapping, by = "PatientID") %>%
  filter(dxBXdays < dxrpdays & bxpos == 1)
bx_1 <- bx_1[,c(8, (2:7))]
bx_1 <- bx_1[order(bx_1$ID),]
  

## psa_1

psa_1 <- psa.ucsf %>% select(-psalessthan) %>%
  left_join(ID_mapping, by = "PatientID") %>%
  filter(dxPSAdays < dxrpdays)
psa_1 <- psa_1[,c(4,2,3)]
psa_1 <- psa_1[order(psa_1$ID),]


## mri_1

mri_1 <- mri.ucsf %>%
  left_join(ID_mapping, by = "PatientID") %>%
  filter(dxMRIdays < dxrpdays)
mri_1 <- mri_1[, c(7,(2:4))]
mri_1 <- mri_1[order(mri_1$ID),]


# Generate dx_2,3, bx_2,3, bxmri_2,3, psa_2,3: add noise to UCSF data

## dx_2
dx_2 <- dx.ucsf
set.seed(011824)
sample_map_2 <- sample(1:nrow(dx_2), 1500, replace = TRUE)
dx_2 <- dx_2[sample_map_2, ]
dx_2$ID <- sample(1:nrow(dx_2), nrow(dx_2))
# Add some noise 
dx_2$age_diag <- round(dx_2$age_diag + rnorm(1500))
dx_2$dxrpdays <- round(dx_2$dxrpdays + rnorm(1500,0,10))
ID_mapping_2 <- dx_2 %>% select(PatientID, ID, dxrpdays)
dx_2 <- dx_2[,c(8, (2:7))]
dx_2 <- dx_2[order(dx_2$ID),]

dx_2$vol <- rep(0, nrow(dx_2))
for (i in 1:nrow(dx_2)){
  dx_2$vol[i] <- rnorm(1, 50 + dx_2$rp[i], 10)
}
dx_2[1300:1500, 2:4] <- NA

## bx_2

bx_2 <- bx.ucsf
bx_2 <- bx_2 %>% 
  left_join(ID_mapping_2, by = "PatientID") %>%
  filter(dxBXdays < dxrpdays & bxpos == 1)
bx_2$gprim <- round(bx_2$gprim + rnorm(nrow(bx_2),0,0.5))
bx_2$gprim[bx_2$gprim > 4] <- 4
bx_2$gprim[bx_2$gprim < 3] <- 3
bx_2$gsecond <- round(bx_2$gsecond + rnorm(nrow(bx_2),0,0.5))
bx_2$gsecond[bx_2$gsecond > 4] <- 4
bx_2$gsecond[bx_2$gsecond < 3] <- 3
bx_2 <- bx_2[,c(8, (2:7))]
bx_2 <- bx_2[order(bx_2$ID),]

## psa_2

psa_2 <- psa.ucsf %>% select(-psalessthan) %>%
  left_join(ID_mapping_2, by = "PatientID") %>%
  filter(dxPSAdays < dxrpdays)

psa_2$psa <- psa_2$psa + rnorm(nrow(psa_2),0,0.5)
psa_2$dxPSAdays <- round(psa_2$dxPSAdays + rnorm(nrow(psa_2),0,10))

psa_2 <- psa_2[,c(4,2,3)]
psa_2 <- psa_2[order(psa_2$ID),]

## mri_2

mri_2 <- mri.ucsf %>%
  left_join(ID_mapping_2, by = "PatientID") %>%
  filter(dxMRIdays < dxrpdays)

mri_2$mriprosvol[!is.na(mri_2$mriprosvol)] <- mri_2$mriprosvol[!is.na(mri_2$mriprosvol)] + 
  rnorm(sum(!is.na(mri_2$mriprosvol)), 0, 5)

mri_2$mripirads <- round(mri_2$mripirads + rnorm(nrow(mri_2),0,0.5))
mri_2$mripirads[mri_2$mripirads > 5] <- 5

mri_2 <- mri_2[, c(7,(2:4))]
mri_2 <- mri_2[order(mri_2$ID),]


## dx_3
dx_3 <- dx.ucsf
set.seed(011824)
sample_map_3 <- sample(1:nrow(dx_3), 1700, replace = TRUE)
dx_3 <- dx_3[sample_map_3, ]
dx_3$ID <- sample(1:nrow(dx_3), nrow(dx_3))
# Add some noise 
dx_3$age_diag <- round(dx_3$age_diag + rnorm(1700))
dx_3$dxrpdays <- round(dx_3$dxrpdays + rnorm(1700,0,10))
ID_mapping_3 <- dx_3 %>% select(PatientID, ID, dxrpdays)
dx_3 <- dx_3[,c(8, (2:7))]
dx_3 <- dx_3[order(dx_3$ID),]
dx_3$vol <- rep(0, nrow(dx_3))
for (i in 1:nrow(dx_3)){
  dx_3$vol[i] <- rnorm(1, 50 + dx_3$rp[i], 10)
}
dx_3[1400:1700, 2:4] <- NA

## bx_3
bx_3 <- bx.ucsf
bx_3 <- bx_3 %>% 
  left_join(ID_mapping_3, by = "PatientID") %>%
  filter(dxBXdays < dxrpdays & bxpos == 1)
bx_3$gprim <- round(bx_3$gprim + rnorm(nrow(bx_3),0,0.5))
bx_3$gprim[bx_3$gprim > 4] <- 4
bx_3$gprim[bx_3$gprim < 3] <- 3
bx_3$gsecond <- round(bx_3$gsecond + rnorm(nrow(bx_3),0,0.5))
bx_3$gsecond[bx_3$gsecond > 4] <- 4
bx_3$gsecond[bx_3$gsecond < 3] <- 3
bx_3 <- bx_3[,c(8, (2:7))]
bx_3 <- bx_3[order(bx_3$ID),]

## psa_3

psa_3 <- psa.ucsf %>% select(-psalessthan) %>%
  left_join(ID_mapping_3, by = "PatientID") %>%
  filter(dxPSAdays < dxrpdays)

psa_3$psa <- psa_3$psa + rnorm(nrow(psa_3),0,0.5)
psa_3$dxPSAdays <- round(psa_3$dxPSAdays + rnorm(nrow(psa_3),0,10))

psa_3 <- psa_3[,c(4,2,3)]
psa_3 <- psa_3[order(psa_3$ID),]

## mri_3

mri_3 <- mri.ucsf %>%
  left_join(ID_mapping_3, by = "PatientID") %>%
  filter(dxMRIdays < dxrpdays)

mri_3$mriprosvol[!is.na(mri_3$mriprosvol)] <- mri_3$mriprosvol[!is.na(mri_3$mriprosvol)] + 
  rnorm(sum(!is.na(mri_3$mriprosvol)), 0, 5)

mri_3$mripirads <- round(mri_3$mripirads + rnorm(nrow(mri_3),0,0.5))
mri_3$mripirads[mri_3$mripirads > 5] <- 5

mri_3 <- mri_3[, c(7,(2:4))]
mri_3 <- mri_3[order(mri_3$ID),]


# Save three synthetic dataset
setwd("/Users/user/Documents/GitHub/prostate-active-surveillance-vDaan/data")
write.csv(dx_1, "dx_1.csv", row.names=FALSE)
write.csv(bx_1, "bx_1.csv", row.names=FALSE)
write.csv(psa_1, "psa_1.csv", row.names=FALSE)
write.csv(mri_1, "mri_1.csv", row.names=FALSE)

write.csv(dx_2, "dx_2.csv", row.names=FALSE)
write.csv(bx_2, "bx_2.csv", row.names=FALSE)
write.csv(psa_2, "psa_2.csv", row.names=FALSE)
write.csv(mri_2, "mri_2.csv", row.names=FALSE)

write.csv(dx_3, "dx_3.csv", row.names=FALSE)
write.csv(bx_3, "bx_3.csv", row.names=FALSE)
write.csv(psa_3, "psa_3.csv", row.names=FALSE)
write.csv(mri_3, "mri_3.csv", row.names=FALSE)