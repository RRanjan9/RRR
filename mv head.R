# Data Entry
head_length_first_son <- c(191, 195, 181, 176, 208, 189, 188, 192, 179, 183, 190, 188, 163, 186, 181, 192)
head_breadth_first_son <- c(155, 149, 148, 144, 157, 150, 152, 152, 158, 147, 159, 151, 137, 153, 140, 154)
head_length_second_son <- c(179, 201, 185, 171, 192, 190, 197, 186, 187, 174, 195, 187, 161, 173, 182, 185)
head_breadth_second_son <- c(145, 152, 149, 142, 152, 149, 159, 151, 148, 147, 157, 158, 130, 148, 146, 152)

# Creating Matrix
data <- cbind(
  x1 = head_length_first_son,
  x2 = head_breadth_first_son,
  x3 = head_length_second_son,
  x4 = head_breadth_second_son
)

# 1. Mean Vector (MLE of μ)
n <- nrow(data)
mean_vector <- apply(data, 2, function(col) sum(col) / n) # Manual mean calculation
cat("Mean Vector (MLE of μ):\n")
print(mean_vector)

# 2. Covariance Matrix (MLE of Σ)
cov_matrix <- matrix(0, ncol = 4, nrow = 4)
for (i in 1:4) {
  for (j in 1:4) {
    cov_matrix[i, j] <- sum((data[, i] - mean_vector[i]) * (data[, j] - mean_vector[j])) / (n - 1)
  }
}
cat("Covariance Matrix (MLE of Σ):\n")
print(cov_matrix)

# 3. Correlation Matrix (ρ)
cor_matrix <- matrix(0, ncol = 4, nrow = 4)
for (i in 1:4) {
  for (j in 1:4) {
    cor_matrix[i, j] <- cov_matrix[i, j] / sqrt(cov_matrix[i, i] * cov_matrix[j, j])
  }
}
cat("Correlation Matrix (ρ):\n")
print(cor_matrix)

# 4. Conditional Distribution Parameters
S11 <- cov_matrix[1:2, 1:2]      # Sub-matrix for (x1, x2)
S12 <- cov_matrix[1:2, 3:4]      # Sub-matrix between (x1, x2) and (x3, x4)
S21 <- t(S12)                    # Transpose of S12
S22 <- cov_matrix[3:4, 3:4]      # Sub-matrix for (x3, x4)
# Solving for S11 Inverse using Manual Inversion (2x2 matrix)
inv_S11 <- matrix(0, 2, 2)
det_S11 <- S11[1, 1] * S11[2, 2] - S11[1, 2] * S11[2, 1]
inv_S11[1, 1] <- S11[2, 2] / det_S11
inv_S11[2, 2] <- S11[1, 1] / det_S11
inv_S11[1, 2] <- -S11[1, 2] / det_S11
inv_S11[2, 1] <- -S11[2, 1] / det_S11

# Conditional Covariance: S22.1 = S22 - S21 * inv(S11) * S12
S22_1 <- S22 - S21 %*% inv_S11 %*% S12
cat("Conditional Covariance Matrix (S22.1):\n")
print(S22_1)

# 5. Partial Correlation Between x3 and x4 Given x1, x2
inv_cov <- solve(cov_matrix) # Manually inverted covariance matrix using built-in solve()
r34.12 <- -inv_cov[3, 4] / sqrt(inv_cov[3, 3] * inv_cov[4, 4])
cat("Partial Correlation r34.12:\n", r34.12, "\n")

# 6. Fisher's Z for Confidence Interval
z_value <- 0.5 * log((1 + r34.12) / (1 - r34.12)) # Fisher's Z-transform
se <- 1 / sqrt(n - 4)
lower <- tanh(z_value - 1.96 * se)
upper <- tanh(z_value + 1.96 * se)
cat("95% Confidence Interval for r34.12:\n")
cat("Lower Bound:", lower, "\nUpper Bound:", upper, "\n")

# 7. Sample Multiple Correlation Coefficient
# Manual Calculation for x3 ~ (x1, x2)
X <- as.matrix(cbind(1, data[, 1:2]))  # Adding intercept term
Y_x3 <- data[, 3]
beta_x3 <- solve(t(X) %*% X) %*% t(X) %*% Y_x3  # Regression Coefficients
Y_hat_x3 <- X %*% beta_x3
RSS_x3 <- sum((Y_x3 - Y_hat_x3)^2)
TSS_x3 <- sum((Y_x3 - mean(Y_x3))^2)
R_x3 <- sqrt(1 - RSS_x3 / TSS_x3)
cat("Multiple Correlation Coefficient (x3 ~ x1, x2):\n", R_x3, "\n")

#Calculation for x4 ~ (x1, x2)
Y_x4 <- data[, 4]
beta_x4 <- solve(t(X) %*% X) %*% t(X) %*% Y_x4
Y_hat_x4 <- X %*% beta_x4
RSS_x4 <- sum((Y_x4 - Y_hat_x4)^2)
TSS_x4 <- sum((Y_x4 - mean(Y_x4))^2)
R_x4 <- sqrt(1 - RSS_x4 / TSS_x4)
cat("Multiple Correlation Coefficient (x4 ~ x1, x2):\n", R_x4, "\n")
