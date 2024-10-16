
library(pROC)
library(ggplot2)
library(gridExtra)

p_eta_est <- apply(out$sims.list$p_eta, c(2,3), mean)
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
    title = sprintf("ROC Curve for (AUC = %.2f)", 0.69), 
    x = "1 - Specificity", 
    y = "Sensitivity"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  annotate("text", x = 0.5, y = 0.1, label = sprintf("AUC: %.2f", 0.69), color = "black", size = 6, fontface = "italic") +
  scale_color_manual(values = c("#00BFC4", "#F8766D"))  # Custom colors

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
    title = sprintf("ROC Curve for (AUC = %.2f)", 0.73), 
    x = "1 - Specificity", 
    y = "Sensitivity"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  annotate("text", x = 0.5, y = 0.1, label = sprintf("AUC: %.2f", 0.73), color = "black", size = 6, fontface = "italic") +
  scale_color_manual(values = c("#00BFC4", "#F8766D"))  # Custom colors

grid.arrange(p1, p2, p3, nrow = 1)


# Diagnostic Misc

#par(mfrow = c(2,2))
#for (i in 1:4){
#  hist(p_eta_est[,i], xlim=c(0,1), main = paste0("P(Cancer State = ",i,")"), xlab="")
#  abline(v = phat[i], col = "red", lwd = 2)
#}

#hist(apply(out$sims.list$eta,2,mean))
#hist(true_rp)
