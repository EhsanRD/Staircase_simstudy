#prepare an example for phd meeting
nsim=1000
S=10
K=10
m=10
ICC=0.2
CAC=0.8
theta=0.15
pow(VarSCcat(S,K,m, ICC, CAC),theta)

set.seed(108159)
res_fit_mat <- replicate(nsim, fitmodels(S, K, m, ICC, CAC, theta))
res_fit_matx_t <-  t(res_fit_mat)
#create a data frame convert res matrix to a data frame
res_fit <- as.data.frame(res_fit_matx_t)
colnames(res_fit) <- c("est_trt_HH", "se_trt_HH", "est_trt_BE", "se_trt_BE", 
                       "est_ICC_HH", "est_ICC_BE", "est_CAC_BE", 
                       "adj_se_KR_HH","adj_ddf_KR_HH","adj_se_KR_NSing_HH", "adj_ddf_KR_NSing_HH",
                       "adj_se_KR_BE","adj_ddf_KR_BE","adj_se_KR_NSing_BE", "adj_ddf_KR_NSing_BE",
                       "adj_se_Sat_HH","adj_ddf_Sat_HH","adj_se_Sat_BE", "adj_ddf_Sat_BE",
                       "w_conv_HH","w_conv_BE","w_other_HH","w_other_BE",
                       "msg_sing_HH","msg_sing_BE","IsSing_HH","IsSing_BE",
                       "err_HH","err_BE")

# Open a PDF file
pdf("histograms_with_true_values_Ex1.pdf", width = 10, height = 8) # Adjust width and height as needed

# List of estimates and true values
est_lst <- list(res_fit$est_trt_HH, res_fit$est_trt_BE, res_fit$est_ICC_HH, 
                res_fit$est_ICC_BE, res_fit$est_CAC_BE)
true_values <- c(theta, theta, ICC, ICC, CAC)  # Replace with actual true values

# List of x-labels for each estimator
x_labels <- c("Estimates for Treatment effect (exchangeable)", 
              "Estimates for Treatment effect (block-exchangeable)", 
              "Estimates for ICC (exchangeable)", "Estimates for ICC (block-exchangeable)",
              "Estimates for CAC (block-exchangeable)")

# Increase top margin
par(mar = c(5, 4, 2, 2))  # Adjust the third element (top margin)

reg <- list(c(0,0.5,0.5,1),c(0.5,1,0.5,1),
            c(0,0.3,0,0.5), c(0.33,0.66,0,0.5),c(0.67,1,0,0.5))

for (i in 1:5) {
  #warning is okay as there is n plot after fifth!
  par(fig= reg[[i]], new=TRUE)
  # Plot histogram
  hist(est_lst[[i]], main = "", 
       xlab = x_labels[i], ylab = "Frequency",cex.axis = 0.8,cex.lab=0.8)
  
  # Add true value to the histogram
  abline(v = true_values[i], col = "red", lwd = 1)
  
  # Calculate new y-coordinate for text
  new_y <- par("usr")[4] * 0.95  # Adjusted to move the text further from the top
  
  # Add text above the vertical line
  #text(true_values[i], new_y, paste("True value:", true_values[i]), col = "black", pos = 1, cex = 0.7)
}

mtext(paste("nsim:",nsim,",","S:",S,",","K:",K,",","m:",m,",",
              "ICC:", ICC,",","CAC:", CAC,",","Theta:", theta),
      side = 3, line =-1.5,outer = TRUE, cex = 0.8)
# Close the PDF file
dev.off()

