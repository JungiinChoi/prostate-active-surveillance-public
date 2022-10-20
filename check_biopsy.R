library(readr)
library(dplyr)
Biopsy_615 <- read_csv("data/Biopsy_6.15.csv")
MRI_615 <- read_csv("data/MRI_6.15.csv")
dim(Biopsy_615[Biopsy_615$Which_Biopsy %in% c("C", "F"),])

#For raw data
length(unique(Biopsy_615$clinical_PTnum))
length(unique(Biopsy_615[Biopsy_615$Which_Biopsy == "D",]$clinical_PTnum))
length(unique(Biopsy_615[Biopsy_615$Which_Biopsy == "C",]$clinical_PTnum))
length(unique(Biopsy_615[Biopsy_615$Which_Biopsy == "F",]$clinical_PTnum))
dt_biopsy <- Biopsy_615 %>% group_by(clinical_PTnum) %>% summarize(n_D = sum(Which_Biopsy == "D"),
                                                                   n_C = sum(Which_Biopsy == "C"),
                                                                   n_F = sum(Which_Biopsy == "F"))
dt_biopsy$clinical_PTnum[dt_biopsy$n_D == 0] #891， removed by Yates
sum(is.na(Biopsy_615$Which_Biopsy)) #941, last Which_Biopsy is NA
Biopsy_615$clinical_PTnum[is.na(Biopsy_615$Which_Biopsy)]
id_hasC <- unique(Biopsy_615$clinical_PTnum[Biopsy_615$Which_Biopsy == "C"])
dt_biopsy_C <- Biopsy_615 %>% 
  filter(clinical_PTnum %in% id_hasC) %>%  
  mutate(Biopsy_Date_f = as.Date(Biopsy_Date, "%m/%d/%y")) %>% 
  arrange(clinical_PTnum, Biopsy_Date_f)

out <- list()
for(i in 1:length(unique(dt_biopsy_C$clinical_PTnum))){
  tmp_dt <- dt_biopsy_C %>% filter(clinical_PTnum == unique(dt_biopsy_C$clinical_PTnum)[i])
  row_C<-which(tmp_dt$Which_Biopsy=="C")
  biop_before_C <- tmp_dt$Which_Biopsy[row_C-1]
  biop_after_C <- tmp_dt$Which_Biopsy[row_C+1]
  out_dt <- c(clinical_PTnum = unique(dt_biopsy_C$clinical_PTnum)[i], 
              biop_before_C = biop_before_C,
              biop_after_C = biop_after_C)
  out[[i]] <- out_dt
}

table(do.call(c, lapply(out, length))) # 3 patient has 2 confirmatory biopsy 
which(do.call(c, lapply(out, length)) == 5)
View(bx.data %>% 
       filter(clinical_PTnum %in% unique(dt_biopsy_C$clinical_PTnum)[which(do.call(c, lapply(out, length)) == 5)]) %>% 
       dplyr::select(clinical_PTnum, Which_Biopsy,mritargetedbx,Biopsy_Date, pgg))
out_dt <- do.call(rbind, out[-which(do.call(c, lapply(out, length)) == 5)])
table(out_dt[,2])
table(out_dt[,3])




# match mri & mri-guided biopsy
load("~/Downloads/prostate-active-surveillance-vDaan/bx.data.RData")
library(zoo)
mri.data <- MRI_615
mri.data$mri_date_f = as.Date(MRI_615$mri_date, "%d%b%Y")
bx.data$biopsy_date_f <- as.Date(bx.data$Biopsy_Date, "%m/%d/%y")
sum(!(unique(mri.data$clinical_ptnum) %in% unique(bx.data$clinical_PTnum))) #100 patients have mri but not in biopsy data
sum(!(unique(mri.data$clinical_ptnum) %in% unique(pt.data$clinical_PTnum))) #100 patients have mri but not in pt data
out_save<-list()
i=9 #this patient has no mritargeted biopsy but does have mri info around biopsy
for(i in 1:length(unique(bx.data$clinical_PTnum))){
#for(i in 10:length(unique(bx.data$clinical_PTnum))){
  print(i)
  bx_dt <- bx.data %>% 
    filter(clinical_PTnum == unique(bx.data$clinical_PTnum)[i]) %>% 
    filter(mritargetedbx==1)
  mri_dt <- mri.data %>%
    filter(clinical_ptnum == unique(bx.data$clinical_PTnum)[i]) %>% 
    arrange(mri_date_f)
  out_match <- bx_dt
  if(nrow(mri_dt) == 0 | nrow(bx_dt) == 0){
   out_match$pirads <- NULL
   out_match$mri_date_f <- NULL
  }else{
    mridate_add <- NULL
    pirads_add <- NULL
    for(j in 1:length(bx_dt$biopsy_date_f)){
      dist_j <- mri_dt$mri_date_f - bx_dt$biopsy_date_f[j]
      #inx_j <- which(abs(dist_j) == min(abs(dist_j)) & dist_j <= 0)
      inx_jsub <- which(dist_j <= 0)
      dist_j_sub <- dist_j[inx_jsub]
      pirads_tmp_sub <- mri_dt$pirads[inx_jsub]
      mridate_tmp <- mri_dt$mri_date_f[inx_jsub]
      
      inx_j <- which(abs(dist_j_sub) == min(abs(dist_j_sub)))
      pirads_tmp <- max(pirads_tmp_sub[inx_j], na.rm = T)
      mridate_tmp <- unique(mridate_tmp[inx_j])
      
      mridate_add <- c(mridate_add, mridate_tmp)
      pirads_add <- c(pirads_add, pirads_tmp)
    }
   out_match$mri_date_f <- mridate_add
   out_match$pirads <- pirads_add
  }
  out_save[[i]] <-  out_match
}
out_nrow <- do.call(c, lapply(out_save, nrow))
out_save_dt <- do.call(rbind, out_save[which(out_nrow != 0)])
bx.data.withmri <- bx.data %>% left_join(out_save_dt)
bx.data.withmri$mri_date_f <- as.Date(bx.data.withmri$mri_date_f)
bx.data.withmri$diff_biopsy_mri <- bx.data.withmri$mri_date_f - bx.data.withmri$biopsy_date_f
summary(abs(as.numeric(bx.data.withmri$diff_biopsy_mri)))
bx.data.withmri %>% filter((as.numeric(diff_biopsy_mri)) < -1000)
hist((as.numeric(bx.data.withmri$diff_biopsy_mri)))


# 3 patients has 2 confirmatory biopsy 
# patient 9 has no mritargeted biopsy but does have mri info around biopsy (check and fix)
# big lag in some biopsy with mri records

#email christian mufaddal and tricia
