#Yates Coley
#rycoley@gmail.com
#2018 February 3, edited 2018 Feb 18, edited 2018 march 16
#This script will load data files and perform initial data checking, tidying, and shaping

#### Libraries required: NONE

### HASHMARK EXPLANATION
# comment explaining code
## workflow description
### header for section of code
#### comment for user (Daan)
##### comment for consideration at a later

### WORKFLOW
## 1. For Daan, specify directory and data file names
## 2. Load libraries 
## 3. Define data checking function. 
## 4. Read in all data files and check unique identifiers 
## 5. Define patient level dataframe (pt.data)
## 6. Create PSA dataframe (psa.data)
## 7. Create biopsy dataframe (bx.data)
## 8. Add treatment data to patient dataframe
## 9. Standardize diagnostic age, prostate volume
## 10. Order patients based on observed true state
## 11. Save datasets.

### 1. For Daan, specify directory and data file names

#### primary directory should contain 3 folders: data, R-scripts, generated-files. 
#This script will pull data from the "data" folder and put an RData database in the "generated-files" folder.
#The complete code will ultimately source several R scripts (all stored in "R-scripts") and save model estimation results in "generated-files" 
# base.location               <- "~/coley/data-cleaning/"
# location.of.data            <- paste0(base.location, "data")
# location.of.r.scripts       <- paste0(base.location, "R-scripts")
# location.of.generated.files <- paste0(base.location, "generated-files")

#### See read.csv() commands below. I have assumed these are .csv files, but this code can be ammended to read in several types of files. 

### 2. Load libraries
### Daan, I didn't want to include an automatic install command, in case you have preferences or limitations in how packages can be installed.
# NO libraries needed at this time

### 3. Define data check function
#this function will terminate program if pathological error identified
data.check <- function(condition, message){
	if(condition==FALSE){print(paste(message, "Program terminated.", sep=" "))}
	stopifnot(condition)}

## 4. Read in all data files and check unique identifiers 
### Data is in a csv-like format, but requires some editting before being 
### workable.
source(paste0(location.of.r.scripts, "/load_data.R"))

### Load csv-files which contain the variable names
dx.names <- read.csv(paste0(location.of.data, "/vars_baseline.csv"), sep = "|")
fu.names <- read.csv(paste0(location.of.data, "/vars_fu.csv"), sep = "|")
tx.names <- read.csv(paste0(location.of.data, "/vars_end.csv"), sep = "|")

dx.data        <- dta_base[, as.character(dx.names$name_gap3)]
names(dx.data) <- dx.names$name_model

fup.data        <- dta_fu[, as.character(fu.names$name_gap3)]
names(fup.data) <- fu.names$name_model

tx.data        <- dta_end[, as.character(tx.names$name_gap3)]
names(tx.data) <- tx.names$name_model

# Remove any rows with missing patient ids
if(sum(is.na(dx.data$id))>0){
  warning(paste0(sum(is.na(dx.data$id)),
                 " rows in the diagnostic dataset were removed because patient id was missing."))}
dx.data <- dx.data[!is.na(dx.data$id),]

if(sum(is.na(fup.data$id))>0){
  warning(paste0(sum(is.na(fup.data$id)),
                 " rows in the follow-up dataset were removed because patient id was missing."))}
fup.data <- fup.data[!is.na(fup.data$id),]

if(sum(is.na(tx.data$id))>0){
  warning(paste0(sum(is.na(tx.data$id)),
                 " rows in the end of AS dataset were removed because patient id was missing."))}
tx.data <- tx.data[!is.na(tx.data$id),]


# Terminate program if there are duplicated unique ids in the patient-level dx data frame
n <- nrow(dx.data)
data.check(condition=as.logical(n==length(unique(dx.data$id))), 
           message="Duplicated PTnum in patient-level diagnosis data frame.")

# Remove any rows without cohort data
if(sum(is.na(dx.data$center))>0){
  warning(paste0(sum(is.na(dx.data$center)), 
                 " individuals in the diagnostic dataset were removed because center was missing."))}

##### Other checks needed?

#number of patients in diagnostic data 
(n<-dim(dx.data)[1])

## 5. Define patient level dataframe (pt.data)

# Just redefined dx.data and pt.data because I will be adding non-diagnostic information to this dataframe
pt.data <- dx.data

# Remove patients who do not meet inclusion criteria
if(sum(pt.data$gs_dx>6 & !is.na(pt.data$gs_dx))>0){
  warning(paste0(sum(pt.data$gs_dx>6 & !is.na(pt.data$gs_dx)),
                 " individuals will not be included in the analysis because diagnostic GS is above 6."))}

if(sum(is.na(pt.data$gs_dx))>0){
  warning(paste0(sum(is.na(pt.data$gs_dx)),
                 " individuals will not be included in the analysis because diagnostic GS is unknown."))}

pt.data <- pt.data[!is.na(pt.data$gs_dx) & pt.data$gs_dx<=6,]
(n<-dim(pt.data)[1])

# Make center variable a factor, name cohort
pt.data$cohort <- as.factor(as.numeric(as.factor(pt.data$center)))
#xtabs(!pt.data$cohort + pt.data$center)
table(pt.data$cohort)

# Define age at diagnosis, remove those without diagnosis age
#pt.data$dx.age <- pt.data$yr_dx - pt.data$yr_birth
names(pt.data)[2] <- "dx.age"
summary(pt.data$dx.age)

if(sum(is.na(pt.data$dx.age))>0){
  warning(paste0(sum(is.na(pt.data$dx.age)),
                 " individuals will not be included in the analysis because diagnostic age is unknown."))}

pt.data <- pt.data[!is.na(pt.data$dx.age),]
(n<-dim(pt.data)[1])

# Remove patients from follow-up and treatment datasets who are excluded from analysis
fup.data <- fup.data[fup.data$id %in% pt.data$id,]
tx.data <- tx.data[tx.data$id %in% pt.data$id,]



### 6. Create PSA dataframe (psa.data)
#get follow-up PSA
psa.data <- fup.data[, c("id", "psa", "days_since_dx")]
names(psa.data) <- c("id", "psa", "days")
psa.data <- psa.data[!is.na(psa.data$psa) & !is.na(psa.data$days),]
(dim(psa.data)[1])

#get PSA observations from one year prior to dx
psa.dx.data <- dx.data[, c("id", "psa_dx")]
names(psa.dx.data) <- c("id", "psa")
psa.dx.data$days <- 0
sum(is.na(psa.dx.data$psa))

#removed those with missing PSA at diagnosis
psa.dx.data <- psa.dx.data[!is.na(psa.dx.data$psa)]

#combine dx and post-dx PSA
psa.data <- rbind(psa.dx.data, psa.data)

#quantify the number of psa observations per patient, remove those with fewer than 2 PSA obs
pt.data$num.psa<-rep(0,n)
for(i in 1:n){pt.data$num.psa[i]<-sum(psa.data$id==pt.data$id[i])}
table(pt.data$num.psa)

if(sum(pt.data$num.psa<2)>0){
  warning(paste0(sum(pt.data$num.psa<2),
                 " individuals will not be included in the analysis because they have fewer than 2 PSA observations."))}

pt.data<-pt.data[pt.data$num.psa>=2,]
(n<-dim(pt.data)[1])

psa.data<-psa.data[psa.data$id %in% pt.data$id,]
length(unique(psa.data$id)) #should be the same as n

#log-transform PSA
sum(psa.data$psa==0)
psa.data$log.psa <- log(psa.data$psa + 0.1)
summary(psa.data$log.psa)


### 7. Create biopsy dataframe (bx.data)
#get follow-up biopsies
bx.data <- fup.data[, c("id", "days_since_dx","num_cores_sampled", 
                        "num_pos_cores", "gleason_sum",
                        "gleason_primary", "gleason_secondary")]
names(bx.data) <- c("id", "days", "num_cores_sampled", "num_pos_cores",
                    "gleason_sum", "gleason_primary", "gleason_secondary")
bx.data <- bx.data[!is.na(bx.data$days),]

bx.data$rm <- rep(0, dim(bx.data)[1])
bx.data$rm[is.na(bx.data$num_pos_cores) & is.na(bx.data$gleason_sum) &
             is.na(bx.data$gleason_primary) & is.na(bx.data$gleason_secondary)] <- 1
bx.data <- bx.data[bx.data$rm==0,]

#assume 12 cores when number of cores sampled is missing
bx.data$num_cores_sampled[is.na(bx.data$num_cores_sampled)] <- 12


## assume gleason primary and secondary = 3 when num_pos_cores=0
if(sum(is.na(bx.data$gleason_primary) & 
       bx.data$num_pos_cores==0 & !is.na(bx.data$num_pos_cores))>0){
  warning(paste0(sum(is.na(bx.data$gleason_primary) & 
                       bx.data$num_pos_cores==0 & !is.na(bx.data$num_pos_cores)),
                 " individuals have a missing primary gleason score but 0 positive cores found. A max-to-date primary Gleason score of 3 is assumed."))}
if(sum(is.na(bx.data$gleason_secondary) & 
       bx.data$num_pos_cores==0 & !is.na(bx.data$num_pos_cores))>0){
  warning(paste0(sum(is.na(bx.data$gleason_secondary) & 
                       bx.data$num_pos_cores==0 & !is.na(bx.data$num_pos_cores)),
          " individuals have a missing secondary gleason score but 0 positive cores found. A max-to-date secondary Gleason score of 3 is assumed."))}

bx.data$gleason_primary[is.na(bx.data$gleason_primary) & 
                          bx.data$num_pos_cores==0 & !is.na(bx.data$num_pos_cores)] <- 3
bx.data$gleason_secondary[is.na(bx.data$gleason_secondary) & 
                          bx.data$num_pos_cores==0 & !is.na(bx.data$num_pos_cores)] <- 3

## assign primary and secondary gleason scores (where missing) based on gleason sum (when available)
## assume gleason primary and secondary=3 when gleason_sum=6
if(sum(is.na(bx.data$gleason_primary) &
       bx.data$gleason_sum==6 & !is.na(bx.data$gleason_sum))>0){
  warning(paste0(sum(is.na(bx.data$gleason_primary) &
                       bx.data$gleason_sum==6 & !is.na(bx.data$gleason_sum)), 
                 " individuals have missing primary Gleason score but Gleason sum is given at 6. Primary Gleason score of 3 is assumed."))}
if(sum(is.na(bx.data$gleason_secondary) &
       bx.data$gleason_sum==6 & !is.na(bx.data$gleason_sum))>0){
  warning(paste0(sum(is.na(bx.data$gleason_secondary) &
                       bx.data$gleason_sum==6 & !is.na(bx.data$gleason_sum)), 
                 " individuals have missing secondary Gleason score but Gleason sum is given at 6. Secondary Gleason score of 3 is assumed."))}

bx.data$gleason_primary[is.na(bx.data$gleason_primary) &
                          bx.data$gleason_sum==6 & !is.na(bx.data$gleason_sum)] <- 3
bx.data$gleason_secondary[is.na(bx.data$gleason_secondary) &
                          bx.data$gleason_sum==6 & !is.na(bx.data$gleason_sum)] <- 3

## assume gleason primary=3 and secondary=4 when gleason_sum=7
#### THIS ASSUMPTION IS CERTAINLY WRONG. We expect about a fifth of biopsy Gleason 7s to be 4+3. 
#### We could do a sensitivity analysis here and assume all are 4+3

#what are the missingness patterns here?
#### IF these show patients with primary missing but not secondary (and vice versa) and Gleason sum of 7, I will write a more complicated code to assign primary and secondary scores
sum(bx.data$gleason_sum==7 & !is.na(bx.data$gleason_sum) & 
      is.na(bx.data$gleason_primary) & !is.na(bx.data$gleason_secondary))
if(sum(bx.data$gleason_sum==7 & !is.na(bx.data$gleason_sum) & 
       is.na(bx.data$gleason_primary) & !is.na(bx.data$gleason_secondary))>0){
  print(bx.data$gleason_secondary[bx.data$gleason_sum==7 & !is.na(bx.data$gleason_sum) & 
                                    is.na(bx.data$gleason_primary) & !is.na(bx.data$gleason_secondary)])}

sum(bx.data$gleason_sum==7 & !is.na(bx.data$gleason_sum) & 
      !is.na(bx.data$gleason_primary) & is.na(bx.data$gleason_secondary))
if(sum(bx.data$gleason_sum==7 & !is.na(bx.data$gleason_sum) & 
       !is.na(bx.data$gleason_primary) & is.na(bx.data$gleason_secondary))>0){
  print(bx.data$gleason_primary[bx.data$gleason_sum==7 & !is.na(bx.data$gleason_sum) & 
                                    !is.na(bx.data$gleason_primary) & is.na(bx.data$gleason_secondary)])}

#assign 3+4, add warning 
if(sum(is.na(bx.data$gleason_primary) &
       bx.data$gleason_sum==7 & !is.na(bx.data$gleason_sum))>0){
  warning(paste0(sum(is.na(bx.data$gleason_primary) &
                       bx.data$gleason_sum==7 & !is.na(bx.data$gleason_sum)), 
                 " individuals have missing primary Gleason score but Gleason sum is given at 7. Primary Gleason score of 3 is assumed."))}
if(sum(is.na(bx.data$gleason_secondary) &
       bx.data$gleason_sum==7 & !is.na(bx.data$gleason_sum))>0){
  warning(paste0(sum(is.na(bx.data$gleason_secondary) &
                       bx.data$gleason_sum==7 & !is.na(bx.data$gleason_sum)), 
                 " individuals have missing secondary Gleason score but Gleason sum is given at 7. Secondary Gleason score of 4 is assumed."))}

bx.data$gleason_primary[is.na(bx.data$gleason_primary) & 
                          bx.data$gleason_sum==7 & !is.na(bx.data$gleason_sum)] <- 3
bx.data$gleason_secondary[is.na(bx.data$gleason_secondary) & 
                          bx.data$gleason_sum==7 & !is.na(bx.data$gleason_sum)] <- 4


## assign 4+4 to patients with missing primary or secondary Gleason score >=8. 4 vs. 5 is inconsequential in this analysis. All are assigned grade group 4+
if(sum(is.na(bx.data$gleason_primary) &
       bx.data$gleason_sum>7 & !is.na(bx.data$gleason_sum))>0){
  warning(paste0(sum(is.na(bx.data$gleason_primary) &
                       bx.data$gleason_sum>7 & !is.na(bx.data$gleason_sum)), 
                 " individuals have missing primary Gleason score but Gleason sum is given greater than 7. Primary Gleason score of 4 is assumed."))}
if(sum(is.na(bx.data$gleason_secondary) &
       bx.data$gleason_sum>7 & !is.na(bx.data$gleason_sum))>0){
  warning(paste0(sum(is.na(bx.data$gleason_secondary) &
                       bx.data$gleason_sum>7 & !is.na(bx.data$gleason_sum)), 
                 " individuals have missing secondary Gleason score but Gleason sum is given greater than 7. Secondary Gleason score of 4 is assumed."))}

bx.data$gleason_primary[is.na(bx.data$gleason_primary) & 
                          bx.data$gleason_sum>7 & !is.na(bx.data$gleason_sum)] <- 4
bx.data$gleason_secondary[is.na(bx.data$gleason_secondary) & 
                          bx.data$gleason_sum>7 & !is.na(bx.data$gleason_sum)] <- 4


# Remove remaining biopsies with missing primary or secondary gleason scores when num_pos_cores > 0 or missing. 
# We don't know what primary and secondary scores are
bx.data$rm[is.na(bx.data$gleason_primary) | is.na(bx.data$gleason_secondary)] <- 1
if(sum(bx.data$rm)>0){
  warning(paste0(sum(bx.data$rm), 
                 " biopsies are missing a primary and/or secondary Gleason score, and one cannot be inferred. These biopsies will be removed from the analysis."))}
bx.data <- bx.data[bx.data$rm==0,]
dim(bx.data)

## Below shouldn't be needed
## Changed the time variable to a factor when recoding the database above. Correct this here:
#bx.data$days <- as.numeric(as.character(bx.data$days))

#quantify the number of bx observations per patient, remove those without any follow-up biopsies
pt.data$num.bx<-rep(0,n)
for(i in 1:n){pt.data$num.bx[i]<-sum(bx.data$id==pt.data$id[i])}
table(pt.data$num.bx)

if(sum(pt.data$num.bx==0)>0){
  warning(paste0(sum(pt.data$num.bx==0),
                 " individuals will not be included in the analysis because they do not have any follow-up biopsies after diagnosis."))}

pt.data<-pt.data[pt.data$num.bx>0,]
(n<-dim(pt.data)[1])

bx.data<-bx.data[bx.data$id %in% pt.data$id,]
length(unique(bx.data$id)) #should be the same as n
(n_bx<-dim(bx.data)[1])

#update psa data with this subset of patients as well
psa.data<-psa.data[psa.data$id %in% pt.data$id,]
(n_psa <- dim(psa.data)[1])
length(unique(psa.data$id)) #should be the same as n

#update tx data
tx.data <- tx.data[tx.data$id %in% pt.data$id,]

#add approximate calendar year of biopsy
bx.data$yr <-vector(length=n_bx)
for(j in 1:n_bx){
  bx.data$yr[j] <- floor(pt.data$yr_dx[pt.data$id==bx.data$id[j]] + bx.data$days[j]/365.25)} # edited
table(bx.data$yr)


#define prognostic grade group of greatest biopsy result to date
bx.data$bx.pgg <- rep(1, n_bx) #default is Gleason 6, grade group 1. will assign this to biopsies with no cancer found as well

### Gleason scores are factors, change to numeric
bx.data$gleason_primary <- as.numeric(as.character(bx.data$gleason_primary))
bx.data$gleason_secondary <- as.numeric(as.character(bx.data$gleason_secondary))

#GS 3+4 = grade group 2
#I use <=3 in case primary or secondary score is 1 or 2
bx.data$bx.pgg[bx.data$gleason_primary<=3 & bx.data$gleason_secondary==4 & 
                 !is.na(bx.data$gleason_primary) & !is.na(bx.data$gleason_secondary)] <- 2
#GS 4+3 = grade group 3
bx.data$bx.pgg[bx.data$gleason_primary==4 & bx.data$gleason_secondary<=3 & 
                 !is.na(bx.data$gleason_primary) & !is.na(bx.data$gleason_secondary)] <- 3
#GS >7 = grade group 4 (or 5)
bx.data$bx.pgg[bx.data$gleason_primary==4 & bx.data$gleason_secondary==4 & 
                 !is.na(bx.data$gleason_primary) & !is.na(bx.data$gleason_secondary)] <- 4
bx.data$bx.pgg[bx.data$gleason_primary>4  & 
                 !is.na(bx.data$gleason_primary) ] <- 4
bx.data$bx.pgg[bx.data$gleason_secondary>4  & 
                 !is.na(bx.data$gleason_secondary) ] <- 4

table(bx.data$bx.pgg)
xtabs(~bx.data$bx.pgg + bx.data$gleason_primary)
xtabs(~bx.data$bx.pgg + bx.data$gleason_secondary)

#binary reclassification
bx.data$rc <- as.numeric(bx.data$bx.pgg>1) # bx.pgg instead of pgg
table(bx.data$rc)

#add final biopsy grade, reclassification, timing to pt.data; remove biopsies that occur after initial reclassification
pt.data$bx.pgg <- pt.data$rc <- pt.data$days.rc <- rep(0,n)
bx.data$rm <- rep(0,n_bx)
for(i in 1:n){
  if(max(bx.data$rc[bx.data$id==pt.data$id[i]])==1){
    pt.data$rc[i]<-1
    pt.data$days.rc[i] <- min(bx.data$days[bx.data$id==pt.data$id[i] & bx.data$rc==1])
    pt.data$bx.pgg[i] <- bx.data$bx.pgg[bx.data$id==pt.data$id[i] & bx.data$days==pt.data$days.rc[i]]
    bx.data$rm[bx.data$id==pt.data$id[i] & bx.data$days>pt.data$days.rc[i]] <- 1 } }
table(pt.data$rc)
table(pt.data$bx.pgg)
xtabs(~pt.data$rc+pt.data$bx.pgg)
xtabs(~pt.data$rc+pt.data$cohort)
xtabs(~pt.data$bx.pgg+pt.data$cohort)

table(bx.data$rm)


##### this may be something we don't want to do in future
bx.data<-bx.data[bx.data$rm==0,]
(n_bx <- dim(bx.data)[1])


### 8. Add treatment data to patient dataframe
#limit treatment data to patients in analysis
tx.data <- tx.data[tx.data$id%in%pt.data$id,]

#limit treatment data to RP
### Currently missing RPs. Suggest to change to grepl("Radical prostatectomy", tx.data$tx_type)
### Following combination also present for instance: Radical prostatectomy and EBRT
tx.data <- tx.data[grepl("Radical prostatectomy", tx.data$tx_type),]
(n_tx <- dim(tx.data)[1])

#look at post-RP gleason grades
table(tx.data$rp_gleason_primary)
table(tx.data$rp_gleason_secondary)
sum(is.na(tx.data$rp_gleason_primary))
sum(is.na(tx.data$rp_gleason_secondary))
sum(is.na(tx.data$rp_gleason_primary) & is.na(tx.data$rp_gleason_secondary))
sum(is.na(tx.data$rp_gleason_primary) & is.na(tx.data$rp_gleason_secondary) 
    & !is.na(tx.data$rp_gleason_sum))

#post-RP grade, recorded in data
tx.data$true.pgg <- rep(NA, n_tx)
tx.data$true.pgg[tx.data$rp_gleason_primary<=3 & tx.data$rp_gleason_secondary<=3 &
                   !is.na(tx.data$rp_gleason_primary) & !is.na(tx.data$rp_gleason_secondary)] <- 1
tx.data$true.pgg[tx.data$rp_gleason_primary<=3 & tx.data$rp_gleason_secondary==4 &
                   !is.na(tx.data$rp_gleason_primary) & !is.na(tx.data$rp_gleason_secondary)] <- 2
tx.data$true.pgg[tx.data$rp_gleason_primary==4 & tx.data$rp_gleason_secondary<=3 &
                   !is.na(tx.data$rp_gleason_primary) & !is.na(tx.data$rp_gleason_secondary)] <- 3
tx.data$true.pgg[tx.data$rp_gleason_primary>4  &
                   !is.na(tx.data$rp_gleason_primary) ] <- 4
tx.data$true.pgg[tx.data$rp_gleason_secondary>4  &
                   !is.na(tx.data$rp_gleason_secondary) ] <- 4
table(tx.data$true.pgg)
sum(is.na(tx.data$true.pgg))

#post-RP grade, when only sum available
#written based on missing Gleason primary. Are any missing primary but not secondary?
sum(is.na(tx.data$rp_gleason_primary) & !is.na(tx.data$rp_gleason_secondary))
sum(!is.na(tx.data$rp_gleason_primary) & is.na(tx.data$rp_gleason_secondary))


## GS sum 6
if(sum(is.na(tx.data$rp_gleason_primary) &
       tx.data$rp_gleason_sum==6 & !is.na(tx.data$rp_gleason_sum))>0){
  warning(paste0(sum(is.na(tx.data$rp_gleason_primary) &
                       tx.data$rp_gleason_sum==6 & !is.na(tx.data$rp_gleason_sum)), 
                 " individuals have missing primary post-RP Gleason score but RP Gleason sum is given as 6. Primary and secondary Gleason score of 3 is assumed."))}
tx.data$true.pgg[tx.data$rp_gleason_sum==6 & !is.na(tx.data$rp_gleason_sum) & 
                   is.na(tx.data$rp_gleason_primary)] <- 1

## GS sum 7, can't determine grade group.
#do we have any?
if(sum(is.na(tx.data$rp_gleason_primary) & 
       tx.data$rp_gleason_sum==7 & !is.na(tx.data$rp_gleason_sum))>0){
  warning(paste0(sum(is.na(tx.data$rp_gleason_primary) & 
                       tx.data$rp_gleason_sum==7 & !is.na(tx.data$rp_gleason_sum)), 
                 " individuals have missing primary post-RP Gleason score and RP Gleason sum==7. Post-RP grade group cannot be determined."))}

## GS sum>7, assign group 4
if(sum(is.na(tx.data$rp_gleason_primary) &
       tx.data$rp_gleason_sum>7 & !is.na(tx.data$rp_gleason_sum))>0){
  warning(paste0(sum(is.na(tx.data$rp_gleason_primary) &
                       tx.data$rp_gleason_sum>7 & !is.na(tx.data$rp_gleason_sum)), 
                 " individuals have missing primary post-RP Gleason score but RP Gleason sum is given as >7. Primary and secondary Gleason score of 4 is assumed."))}
tx.data$true.pgg[tx.data$rp_gleason_sum>7 & !is.na(tx.data$rp_gleason_sum) & 
                   is.na(tx.data$rp_gleason_primary)] <- 4


#add surgery indicator and result to pt.data
pt.data$rp <- rep(0,n)
pt.data$true.pgg <- rep(NA, n)
for(i in 1:n){
  if(pt.data$id[i]%in%tx.data$id){
    pt.data$rp[i] <- 1
    pt.data$true.pgg[i] <- tx.data$true.pgg[tx.data$id==pt.data$id[i]] } }
table(pt.data$rp)
table(pt.data$true.pgg)
sum(pt.data$rp==1 & is.na(pt.data$true.pgg))

xtabs(~pt.data$rp + pt.data$cohort)
xtabs(~pt.data$rp + pt.data$true.pgg)
xtabs(~pt.data$true.pgg + pt.data$cohort)

xtabs(~pt.data$true.pgg + pt.data$bx.pgg)
xtabs(!pt.data$rp + pt.data$bx.pgg)



### 9. Standardize diagnostic age, prostate volume

##standarize age
pt.data$dx.age.std <- scale(pt.data$dx.age, center=T, scale=T)
summary(pt.data$dx.age.std)

#put it in PSA dataframe
psa.data$dx.age.std <- rep(0, n_psa)
for(i in 1:n){
  psa.data$dx.age.std[psa.data$id==pt.data$id[i]] <- pt.data$dx.age.std[i] }


##standardize prostate volume
pt.data$vol.std[!is.na(pt.data$vol)] <- scale(pt.data$vol[!is.na(pt.data$vol)], 
                                              center=mean(pt.data$vol, na.rm=T),
                                              scale=sd(pt.data$vol, na.rm=T))
pt.data$vol.std[is.na(pt.data$vol)] <- 0 #assume average prostate volume when unknown
summary(pt.data$vol.std)

#put in PSA dataframe
psa.data$vol.std <- rep(0, n_psa)
for(i in 1:n){
  psa.data$vol.std[psa.data$id==pt.data$id[i]] <- pt.data$vol.std[i] }

#put in BX dataframe
bx.data$vol.std<-rep(0, n_bx)
for(i in 1:n){
  bx.data$vol.std[bx.data$id==pt.data$id[i]] <- pt.data$vol.std[i] }

### 10. Order patients based on observed true state
#### The need to order patient data in this way is an artifact of JAGS indexing requirements
ordered<-order(pt.data$true.pgg)
pt.data<-pt.data[ordered,]
pt.data$subj<-c(1:n)

#add this unique identifier to PSA and BX data
psa.data$subj<-rep(0,n_psa)
bx.data$subj<-rep(0,n_bx)
for(i in 1:n){
	psa.data$subj[psa.data$id==pt.data$id[i]] <- i
	bx.data$subj[bx.data$id==pt.data$id[i]] <- i }

### 11. Save data
save(pt.data, psa.data, bx.data,  
     file=paste0(location.of.generated.files,"/AS-data-shaping-work-space.RData"))




### 12. Additional descriptives and checks
summary(pt.data)
dim(pt.data)[1]
length(unique(pt.data$subj))
table(pt.data$cohort)



summary(psa.data)
dim(psa.data)[1]
length(unique(psa.data$subj))

summary(bx.data)
dim(bx.data)[1]
table(bx.data$bx.pgg)

