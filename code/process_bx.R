base.location <- "/Users/zitongwang/Downloads/prostate-active-surveillance-vDaan/"
location.of.data <- paste0(base.location, "data")
location.of.r.scripts <- paste0(base.location, "cleaned_code_data2022")
location.of.generated.files <- paste0(base.location, "cleaned_code_data2022/generated-files")


name.of.pt.data <- "Demographics_10.29.22.csv" 
name.of.bx.data <- "Biopsy_11.7.22.csv" 
name.of.psa.data <- "PSA_10.29.22.csv" 
name.of.tx.data <- "Treatments_10.29.22.csv" 
name.of.mri.data <- "MRI_10.29.22.csv" 
name.of.mri.findings.data <- "MRI_findings_10.29.22.csv"
name.of.targeted.pathology.data <- "Targeted_pathology_10.29.22.csv"

target.path.data <- read_csv(paste(location.of.data, name.of.targeted.pathology.data, sep = "/"))
mri.find.data <- read_csv(paste(location.of.data, name.of.mri.findings.data, sep = "/"))

bx.data <- read_csv(paste(location.of.data, name.of.bx.data, sep = "/"))

mri.data <- read_csv(paste(location.of.data, name.of.mri.data, sep = "/"))
mri.problem <- problems()
mri.data[mri.problem$row-1,mri.problem$col] <- as.numeric(gsub(",", "", mri.problem$actual))


mri.data.new <- mri.data %>% left_join(mri.find.data, by = "mri_id")
#assign all missing pirads but no_measurable_disease=1 to have pirads 1
mri.data.new<- mri.data.new %>% mutate(pirads_update  = 
                                         ifelse(is.na(pirads) & no_measurable_disease == 1, 
                                                1,
                                                pirads))

target.path.data.short <- target.path.data %>% group_by(Biopsyid) %>% summarize(max_pirads = max(pirads))
bx.data.new<- bx.data %>% left_join(target.path.data.short, by = "Biopsyid")
nrow(bx.data) == nrow(bx.data.new)

write.csv(bx.data.new, file = "data/Biopsy_processed.csv")


name.of.bx.data <- "Biopsy_processed.csv" 

