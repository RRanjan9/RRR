# Set parameters
shape <- 2
scale <- 5
n <- 100

# Generate 100 observations from Gamma distribution
gamma_data <- rgamma(n, shape = shape, scale = scale)

# Gaussian Kernel function based on the formula
gaussian_kernel <- function(u) {
  return((1 / sqrt(2 * pi)) * exp(-0.5 * u^2))
}

# Function to estimate density using the given formula
kernel_density_estimate <- function(x_values, data, h) {
  n <- length(data)
  density_estimates <- numeric(length(x_values))
  
  for (i in 1:length(x_values)) {
    sum_kernels <- 0
    for (j in 1:n) {
      sum_kernels <- sum_kernels + gaussian_kernel((x_values[i] - data[j]) / h)
    }
    density_estimates[i] <- (1 / (n * h)) * sum_kernels
  }
  
  return(density_estimates)
}

# Define the range of x values for plotting the density
x_values <- seq(min(gamma_data), max(gamma_data), length.out = 100)

# Bandwidths: 0.1, 0.2, 0.5
bw_0.1 <- kernel_density_estimate(x_values, gamma_data, 0.1)
bw_0.2 <- kernel_density_estimate(x_values, gamma_data, 0.2)
bw_0.5 <- kernel_density_estimate(x_values, gamma_data, 0.5)

# Print the estimated density values for each bandwidth
cat("Kernel Density Estimates for bandwidth 0.1:\n")
print(bw_0.1)

cat("\nKernel Density Estimates for bandwidth 0.2:\n")
print(bw_0.2)

cat("\nKernel Density Estimates for bandwidth 0.5:\n")
print(bw_0.5)

# Plot the kernel density estimates
plot(x_values, bw_0.1, type = "l", col = "red", lwd = 2, main = "Kernel Density Estimates (Gamma)", xlab = "x", ylab = "Density")
lines(x_values, bw_0.2, col = "blue", lwd = 2)
lines(x_values, bw_0.5, col = "green", lwd = 2)

# Add a legend
legend("topright", legend = c("bw = 0.1", "bw = 0.2", "bw = 0.5"), col = c("red", "blue", "green"), lwd = 2)