### Define directories, file names
base.location <- workdir
location.of.data <- paste0(base.location, "/data")
location.of.r.scripts <- paste0(base.location, "/code")
location.of.generated.files <- paste0(base.location, "/generated-files")
location.of.generated.folder = paste(location.of.generated.files, "/hier",sep="")

library(pROC)
library(ggplot2)
library(gridExtra)

# Number of Institutions
J <- 3

p_eta_est <- list(length = J)
for(j in 1:J){
  p_eta_est[[j]] <- read.csv(paste(location.of.generated.folder, "/jags-prediction-p_eta_", 
                                   j,"-", mri_role,".csv",sep=""))[,-1]
}

# AUC dataframe for PGG > 1, PGG > 2, PGG > 3
# Plot ROC Curve

dx_list <- list(length = J)
AUC_independent <- matrix(0,nrow = J, ncol = 3)
colnames(AUC_independent) <- c("> 1", "> 2", "> 3")
plots_lst <- list()
ind <- 1

for (i in 1:J){
  dx_data <- read.csv(paste0(location.of.data,"/dx_",i,".csv"))
  dx_list[[i]] <- dx_data
  true_ind <- !is.na(dx_data$rp)
  true_rp <- dx_data$rp[true_ind]
  
  for (k in 1:3){
    true_class <- (true_rp > k)
    predicted_prob <- if(k==3){
      predicted_prob <- p_eta_est[[i]][,4][1:nrow(dx_data)][true_ind]
    }else{rowSums(p_eta_est[[i]][,-(1:k)])[1:nrow(dx_data)][true_ind]}
    
    # Calculate the ROC curve and AUC
    roc_obj <- roc(true_class, predicted_prob)
    
    # Extract data for plotting ROC curve
    roc_data <- data.frame(
      specificity = rev(roc_obj$specificities),
      sensitivity = rev(roc_obj$sensitivities)
    )
    
    
    # Plot ROC curve using ggplot2
    plots_lst[[ind]] <- ggplot(roc_data, aes(x = specificity, y = sensitivity)) +
      geom_line(color = "#00BFC4", linewidth = 1.5) +   # ROC curve
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +  # Diagonal line
      theme_minimal(base_size = 15) +  # Minimal theme with larger text
      labs(
        title = paste0("Center ", i, ": PGG > ",k),
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
    
    AUC_independent[i,k] <- auc(roc_obj)
    ind <- ind + 1
  }
  
}

grid.arrange(grobs = plots_lst, nrow = J)

# Matrix of Independent AUC

AUC_independent

