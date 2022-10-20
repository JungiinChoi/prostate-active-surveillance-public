# check for drop in gleason scores 
# for running data-load-check-and-shaping.R till line 435
library(dplyr)
n.bx.data <- bx.data %>% 
  group_by(clinical_PTnum) %>% 
  filter(!is.na(pgg))%>% 
  summarise(nbiopsy = n(),
            index = 1:nbiopsy) 
n.bx.data <- cbind(n.bx.data, bx.data %>% filter(!is.na(pgg)) %>% dplyr::select(Biopsy_Date,pgg))

max.bx.data <- n.bx.data %>%
  group_by(clinical_PTnum) %>% 
  slice_max(pgg) %>% 
  dplyr::select(clinical_PTnum, index, pgg) %>% 
  summarize(loc_maxbiopsy = last(index))
last.bx.data <- n.bx.data %>% dplyr::select(clinical_PTnum, nbiopsy) %>% unique()
check.bx.data <- last.bx.data %>% full_join(max.bx.data)

cpat <- check.bx.data$clinical_PTnum[which(check.bx.data$loc_maxbiopsy < check.bx.data$nbiopsy)]
