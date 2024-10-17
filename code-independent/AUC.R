
### Define directories, file names
base.location <- workdir
location.of.data <- paste0(base.location, "/data-independent")
location.of.r.scripts <- paste0(base.location, "/code-independent")
location.of.generated.files <- paste0(base.location, "/generated-files")
location.of.generated.folder = paste(location.of.generated.files, "/indep",sep="")

library(pROC)
library(ggplot2)
library(gridExtra)

p_eta_est <- read.csv(paste(location.of.generated.folder, "/jags-prediction-p_eta-", 
                            mri_role,".csv",sep=""))[,-1]
dx_data <- read.csv(paste0(location.of.data,"/dx.csv"))
p_eta_est <- apply(p_eta_hat, c(2,3), mean)
true_ind <- !is.na(dx_data$rp)
true_rp <- dx_data$rp[true_ind]

k <- 1



true_class <- (true_rp > k)
predicted_prob <- (rowSums(p_eta_est[,-(1:k)]))[true_ind]
# Calculate the ROC curve and AUC
roc_obj <- roc(true_class, predicted_prob)

# Extract data for plotting ROC curve
roc_data <- data.frame(
  specificity = rev(roc_obj$specificities),
  sensitivity = rev(roc_obj$sensitivities)
)

# Plot ROC curve using ggplot2
p1 <- ggplot(roc_data, aes(x = specificity, y = sensitivity)) +
  geom_line(color = "#00BFC4", linewidth = 1.5) +   # ROC curve
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +  # Diagonal line
  theme_minimal(base_size = 15) +  # Minimal theme with larger text
  labs(
    title = sprintf("ROC Curve for (AUC = %.2f)", auc(roc_obj)), 
    x = "1 - Specificity", 
    y = "Sensitivity"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  annotate("text", x = 0.5, y = 0.1, label = sprintf("AUC: %.2f", auc(roc_obj)), color = "black", size = 6, fontface = "italic") +
  scale_color_manual(values = c("#00BFC4", "#F8766D"))  # Custom colors
AUC_indepentent <- auc(roc_obj)

k <- 2

true_class <- (true_rp > k)
predicted_prob <- (rowSums(p_eta_est[,-(1:k)]))[true_ind]


# Calculate the ROC curve and AUC
roc_obj <- roc(true_class, predicted_prob)

# Extract data for plotting ROC curve
roc_data <- data.frame(
  specificity = rev(roc_obj$specificities),
  sensitivity = rev(roc_obj$sensitivities)
)

# Plot ROC curve using ggplot2
p2 <- ggplot(roc_data, aes(x = specificity, y = sensitivity)) +
  geom_line(color = "#00BFC4", linewidth = 1.5) +   # ROC curve
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +  # Diagonal line
  theme_minimal(base_size = 15) +  # Minimal theme with larger text
  labs(
    title = sprintf("ROC Curve for (AUC = %.2f)", auc(roc_obj)), 
    x = "1 - Specificity", 
    y = "Sensitivity"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  annotate("text", x = 0.5, y = 0.1, label = sprintf("AUC: %.2f", auc(roc_obj)), color = "black", size = 6, fontface = "italic") +
  scale_color_manual(values = c("#00BFC4", "#F8766D"))  # Custom colors
AUC_indepentent <- c(AUC_indepentent, auc(roc_obj))

k <- 3

true_class <- (true_rp > k)
predicted_prob <- p_eta_est[,4][true_ind]

# Calculate the ROC curve and AUC
roc_obj <- roc(true_class, predicted_prob)

# Extract data for plotting ROC curve
roc_data <- data.frame(
  specificity = rev(roc_obj$specificities),
  sensitivity = rev(roc_obj$sensitivities)
)

# Plot ROC curve using ggplot2
p3 <- ggplot(roc_data, aes(x = specificity, y = sensitivity)) +
  geom_line(color = "#00BFC4", linewidth = 1.5) +   # ROC curve
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +  # Diagonal line
  theme_minimal(base_size = 15) +  # Minimal theme with larger text
  labs(
    title = sprintf("ROC Curve for Independent Model (AUC = %.2f)", auc(roc_obj)), 
    x = "1 - Specificity", 
    y = "Sensitivity"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  annotate("text", x = 0.5, y = 0.1, label = sprintf("AUC: %.2f", auc(roc_obj)), color = "black", size = 6, fontface = "italic") +
  scale_color_manual(values = c("#00BFC4", "#F8766D"))  # Custom colors

AUC_indepentent <- c(AUC_indepentent, auc(roc_obj))

grid.arrange(p1, p2, p3, nrow = 1)

# List of Independent AUC

names(AUC_indepentent) <- c("> 1", "> 2", "> 3")
AUC_indepentent


# Diagnostic Misc

#par(mfrow = c(2,2))
#for (i in 1:4){
#  hist(p_eta_est[,i], xlim=c(0,1), main = paste0("P(Cancer State = ",i,")"), xlab="")
#  abline(v = phat[i], col = "red", lwd = 2)
#}

#hist(apply(out$sims.list$eta,2,mean))
#hist(true_rp)


library(pROC)
library(ggplot2)
library(gridExtra)

p_eta_est <- apply(out$sims.list$p_eta, c(2,3), mean)
true_ind <- !is.na(dx_data$rp)
true_rp <- dx_data$rp[true_ind]
AUC_independent <- matrix(0, nrow = 3, ncol = J)

for (k in 1:3){
  true_class <- (true_rp > k)
  predicted_prob <- if(k==3){
    predicted_prob <- p_eta_est[,4][true_ind]
  }else{rowSums(p_eta_est[,-(1:k)]))[true_ind]}
  
  # Calculate the ROC curve and AUC
  roc_obj <- roc(true_class, predicted_prob)
  
  # Extract data for plotting ROC curve
  roc_data <- data.frame(
    specificity = rev(roc_obj$specificities),
    sensitivity = rev(roc_obj$sensitivities)
  )
  
  plots_lst <- list()
  
  # Plot ROC curve using ggplot2
  plots_lst[[k]] <- ggplot(roc_data, aes(x = specificity, y = sensitivity)) +
    geom_line(color = "#00BFC4", linewidth = 1.5) +   # ROC curve
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +  # Diagonal line
    theme_minimal(base_size = 15) +  # Minimal theme with larger text
    labs(
      title = sprintf("ROC Curve for (AUC = %.2f)", auc(roc_obj)), 
      x = "1 - Specificity", 
      y = "Sensitivity"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12)
    ) +
    annotate("text", x = 0.5, y = 0.1, label = sprintf("AUC: %.2f", auc(roc_obj)), color = "black", size = 6, fontface = "italic") +
    scale_color_manual(values = c("#00BFC4", "#F8766D"))  # Custom colors
  AUC_indepentent <- c(AUC_indepentent, auc(roc_obj))
}

grid.arrange(grobs = plots_lst, nrow = 1)

grid.arrange(p1, p2, p3, nrow = 1)

# List of Independent AUC

names(AUC_indepentent) <- c("> 1", "> 2", "> 3")
AUC_indepentent


# Diagnostic Misc

#par(mfrow = c(2,2))
#for (i in 1:4){
#  hist(p_eta_est[,i], xlim=c(0,1), main = paste0("P(Cancer State = ",i,")"), xlab="")
#  abline(v = phat[i], col = "red", lwd = 2)
#}

#hist(apply(out$sims.list$eta,2,mean))
#hist(true_rp)


