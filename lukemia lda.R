# Given data
lifetimes <- c(22.3, 26.8, 30.3, 31.9, 32.1, 33.3, 33.7, 33.9, 34.7, 36.1, 36.4, 36.5, 
               36.6, 37.1, 37.6, 38.2, 38.5, 38.7, 38.7, 38.9, 38.9, 39.1, 41.1, 41.1, 
               41.4, 42.4, 43.6, 43.8, 44.0, 45.3, 45.8, 50.4, 51.3, 51.4, 51.5)
# Log-transform the data
log_lifetimes <- log(lifetimes)

# Estimate parameters (mu and sigma) using MLE
mu_hat <- mean(log_lifetimes)
sigma_hat <- sd(log_lifetimes)

# Calculate Mean and Median Time to Failure
mean_time_to_failure <- exp(mu_hat + (sigma_hat^2 / 2))
median_time_to_failure <- exp(mu_hat)

# Survival Function
survival_function <- function(t) {
  1 - pnorm((log(t) - mu_hat) / sigma_hat)
}
# Hazard Function
hazard_function <- function(t) {
  dnorm((log(t) - mu_hat) / sigma_hat) / (t * survival_function(t))
}
# Generate a sequence of times for plotting
time_seq <- seq(min(lifetimes), max(lifetimes), length.out = 100)
# Survival and Hazard values
survival_vals <- survival_function(time_seq)
hazard_vals <- hazard_function(time_seq)

# Plot Survival Curve
plot(time_seq, survival_vals, type = "l", col = "blue", lwd = 2,
     xlab = "Time (weeks)", ylab = "Survival Probability",
     main = "Survival Curve")
# Plot Hazard Curve
plot(time_seq, hazard_vals, type = "l", col = "red", lwd = 2,
     xlab = "Time (weeks)", ylab = "Hazard Rate",
     main = "Hazard Curve")
# results
cat("Estimated Parameters:\n")
cat("Mu (mean of log-lifetimes):", mu_hat, "\n")
cat("Sigma (std. dev. of log-lifetimes):", sigma_hat, "\n\n")
cat("Mean Time to Failure (MTTF):", mean_time_to_failure, "weeks\n")
cat("Median Time to Failure:", median_time_to_failure, "weeks\n")
