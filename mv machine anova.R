# Input Data
data <- matrix(
  c(15, 14, 19, 18, 
    17, 12, 20, 16, 
    16, 18, 16, 17, 
    16, 16, 15, 15), 
  nrow = 4, byrow = TRUE
)

rownames(data) <- c("Machinist1", "Machinist2", "Machinist3", "Machinist4")
colnames(data) <- c("Machine1", "Machine2", "Machine3", "Machine4")

# Total Observations
N <- length(data) # Total number of observations
n_row <- nrow(data) # Number of Machinists
n_col <- ncol(data) # Number of Machines

# Step 1: Grand Mean
grand_mean <- mean(data)

# Step 2: Total Sum of Squares (SST)
SST <- sum((data - grand_mean)^2)

# Step 3: Row Sum of Squares (SSR - for Machinists)
row_means <- rowMeans(data)
SSR <- n_col * sum((row_means - grand_mean)^2)

# Step 4: Column Sum of Squares (SSC - for Machines)
col_means <- colMeans(data)
SSC <- n_row * sum((col_means - grand_mean)^2)

# Step 5: Interaction Sum of Squares (SSI)
interaction_component <- 0
for (i in 1:n_row) {
  for (j in 1:n_col) {
    interaction_component <- interaction_component + 
      (data[i, j] - row_means[i] - col_means[j] + grand_mean)^2
  }
}
SSI <- interaction_component

# Step 6: Residual/Error Sum of Squares (SSE)
SSE <- SST - SSR - SSC - SSI

# Step 7: Mean Squares
MSR <- SSR / (n_row - 1)    # Mean Square for Rows (Machinists)
MSC <- SSC / (n_col - 1)    # Mean Square for Columns (Machines)
MSI <- SSI / ((n_row - 1) * (n_col - 1))  # Mean Square for Interaction
MSE <- SSE / ((n_row - 1) * (n_col - 1))  # Mean Square for Error

# Step 8: F-statistics
F_Rows <- MSR / MSE        # F-statistic for Rows
F_Columns <- MSC / MSE     # F-statistic for Columns
F_Interaction <- MSI / MSE # F-statistic for Interaction

# Step 9: Display Results
cat("Two-Way ANOVA Results:\n")
cat("-------------------------------------------------\n")
cat("SST (Total Sum of Squares):", SST, "\n")
cat("SSR (Machinist Sum of Squares):", SSR, "\n")
cat("SSC (Machine Sum of Squares):", SSC, "\n")
cat("SSI (Interaction Sum of Squares):", SSI, "\n")
cat("SSE (Error Sum of Squares):", SSE, "\n\n")

cat("Mean Squares:\n")
cat("MSR (Rows):", MSR, "\n")
cat("MSC (Columns):", MSC, "\n")
cat("MSI (Interaction):", MSI, "\n")
cat("MSE (Error):", MSE, "\n\n")

cat("F-statistics:\n")
cat("F for Rows (Machinists):", F_Rows, "\n")
cat("F for Columns (Machines):", F_Columns, "\n")
cat("F for Interaction:", F_Interaction, "\n")
