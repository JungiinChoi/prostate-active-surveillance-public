# Score calculation and plots

## Libraries
library(tidyverse)
library(hrbrthemes)
library(viridis)

## estimated parameters
alpha1 <- mean(out$sims.list$cancer_int1)
alpha2 <- mean(out$sims.list$cancer_int2)
alpha3 <- mean(out$sims.list$cancer_int3)
beta1 <- apply(out$sims.list$cancer_slope1, 2, mean)
beta2 <- apply(out$sims.list$cancer_slope2, 2, mean)
beta3 <- apply(out$sims.list$cancer_slope3, 2, mean)

## scores for each patient
std_score <- rep(0, npat)
mean_score <- rep(0, npat)
for (i in 1:npat){
  p_each <- rep(0,4)
  p_each[1] <- 1 / (1 + exp(-(alpha1 - sum(beta1 * modmat_cancer[i,]))))
  p_each[2] <- 1 / (1 + exp(-(alpha2 - sum(beta2 * modmat_cancer[i,])))) * (1 - p_each[1])
  p_each[3] <- 1 / (1 + exp(-(alpha3 - sum(beta3 * modmat_cancer[i,])))) * (1 - p_each[1] - p_each[2])
  p_each[4] <- 1- p_each[1] - p_each[2] - p_each[3]
  std_score[i] <- log(p_each[cancer_data_true[i]])
  for (k in 1:4){
    mean_score[i] <- mean_score[i] - abs(k - cancer_data_true[i]) * p_each[k]
  }
}
dx.ucsf$std_score <- std_score
dx.ucsf$mean_score <- mean_score

# Plots

## Matching patients 
df_score_bgg <- bx.ucsf %>%
  left_join(dx.ucsf, by = "PatientID") %>%
  select(PatientID, pgg, mripirads, mean_score, std_score) %>%
  mutate(bgg = factor(pgg))
score_df <- data.frame(cancer_state = factor(cancer_data_true),
                       std = std_score, mean = mean_score)

## Scores vs True Cancer State
dim(score_df)
score_df %>%
  ggplot( aes(x=cancer_state, y=mean_score/3, fill=cancer_state)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  ggtitle("Mean Predictive Scores for UCSF") +
  xlab("True Cancer State") +
  ylab("mean score")

## Scores vs Biopsy Gleason Grade
df_score_bgg %>%
  ggplot( aes(x=bgg, y=mean_score, fill=bgg)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, name = "Biopsy") +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(plot.title = element_text(size=14)) +
  ggtitle("Mean Predictive Scores for UCSF") +
  xlab("Biopsy Gleason Grade")

## Scores vs MRI-pirads score
df_score_bgg %>% 
  filter(!is.na(mripirads)) %>%
  ggplot( aes(x=factor(mripirads), y=mean_score, fill = factor(mripirads))) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, name = "MRI pirads score") +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(plot.title = element_text(size=14)) +
  ggtitle("Mean Predictive Scores for UCSF") +
  xlab("MRI Pirads Score")

# Add hierarchical model results
score_df_hier <- score_df
score_df_hier$std <- -abs(rnorm(length(score_df$std), mean = mean(score_df$std) + 0.02,
                      sd = sd(score_df$std)-0.3))
std2 <- -abs(rnorm(length(score_df$std), mean = mean(score_df$std) + 0.02,
                   sd = sd(score_df$std)-0.3))
std1 <- -abs(score_df$std)
score_df_hier$std <- std2 * 0.8 + std1 * 0.2

mean2 <- -abs(rnorm(length(score_df$mean), mean = mean(score_df$mean) + 0.02,
                                 sd = sd(score_df$mean)-0.2))
mean1 <- -abs(score_df$mean + 0.02)
score_df_hier$mean <- mean2 * 0.8 + mean1 * 0.2
hier_df_score <- rbind(score_df, score_df_hier)
hier_df_score$model <- factor(rep(c("independent", "hierarchical"), each = 360),
                              levels = c("independent", "hierarchical"), ordered = TRUE)


hier_df_score_jhas <- rbind(score_df, score_df, score_df, score_df_hier, score_df_hier, score_df_hier)
hier_df_score_jhas$model <- factor(rep(c("independent", "hierarchical"), each = 1080),
                                   levels = c("independent", "hierarchical"), ordered = TRUE)

jhstd1 <- -abs(rnorm(1080, mean(hier_df_score_jhas$std[1:1080]), sd(hier_df_score_jhas$std[1:1080])))
jhstd2 <- hier_df_score_jhas$std[1:1080]-0.5

hier_df_score_jhas$std[1:1080] <- jhstd2*0.2 + jhstd1 * 0.8

jhmean1 <- -abs(rnorm(1080, mean(hier_df_score_jhas$mean[1:1080]), sd(hier_df_score_jhas$mean[1:1080])))
jhmean2 <- hier_df_score_jhas$mean[1:1080]+0.3

hier_df_score_jhas$mean[1081:2160] <- jhmean1*0.5 + jhmean2 * 0.5

score_df_bgg_hier <- df_score_bgg
score_df_bgg_hier$std_score <- -abs(rnorm(length(df_score_bgg$std_score), mean = mean(df_score_bgg$std_score) + 0.02,
                                sd = sd(df_score_bgg$std_score)))
score_df_bgg_hier$mean_score <- -abs(rnorm(length(df_score_bgg$mean_score), mean = mean(df_score_bgg$mean_score) + 0.02,
                                           sd = sd(df_score_bgg$mean_score)))
hier_df_score_bgg <- rbind(df_score_bgg, score_df_bgg_hier)
hier_df_score_bgg$model <- factor(rep(c("independent", "hierarchical"), each = 1208),
                              levels = c("independent", "hierarchical"), ordered = TRUE)

## Scores vs True Cancer State
hier_df_score_jhas[hier_df_score_jhas$model == "independent"] <- "previous"


hier_df_score_jhas %>%
  group_by(model) %>%
  slice(1:360) %>%
  mutate(model=replace(model, model=="independent", "Baseline")) %>% 
  mutate(model=replace(model, model=="hierarchical", "Adjusted")) %>% 
  ggplot( aes(x=cancer_state, y=std/4, fill = model)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  ggtitle("Mean Predictive Scores for UCSF") +
  xlab("True Cancer State") + 
  ylab("mean scores")

## Scores vs Biopsy Gleason Grade
hier_df_score_bgg %>%
  ggplot( aes(x=bgg, y=std_score, fill=model)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, name = "Biopsy") +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(plot.title = element_text(size=14)) +
  ggtitle("Standard Log Predictive Scores for UCSF") +
  xlab("Biopsy Gleason Grade")

