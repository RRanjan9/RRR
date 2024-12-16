death_times <- c(7.35, 8.69, 8.80, 9.63, 9.63, 9.89, 9.98, 10.24, 10.36, 10.37, 
                 10.48, 11.33, 11.39, 12.02, 13.12)  # Complete data
censored_times <- rep(20, 10)  # Right-censored data
# Combine data
all_times <- c(death_times, censored_times)
status <- c(rep(1, length(death_times)), rep(0, length(censored_times)))  # 1=death, 0=censored

# Negative log-likelihood function
neg_log_likelihood <- function(params) {
  alpha <- params[1]  # Shape parameter
  beta <- params[2]   # Scale parameter
  if (alpha <= 0 || beta <= 0) return(Inf) 
  log_f <- log(alpha / beta) + (alpha - 1) * log(death_times / beta) - (death_times / beta)^alpha
  log_S <- -(censored_times / beta)^alpha
  total_log_likelihood <- sum(log_f) + sum(log_S)
  return(-total_log_likelihood)  # Negative log-likelihood for minimization
}
# Initial guesses for alpha and beta
initial_guess <- c(1, 15)

# MLE using optim
result <- optim(par = initial_guess, fn = neg_log_likelihood, method = "L-BFGS-B", 
                lower = c(1e-5, 1e-5))
mle_alpha <- result$par[1]
mle_beta <- result$par[2]

# Survival function
survival_function <- function(t) {
  exp(-(t / mle_beta)^mle_alpha)
}
# Hazard function
hazard_function <- function(t) {
  (mle_alpha / mle_beta) * (t / mle_beta)^(mle_alpha - 1)
}
# Generate time sequence for plotting
time_seq <- seq(0, 25, length.out = 100)
# Survival and hazard values
survival_vals <- survival_function(time_seq)
hazard_vals <- hazard_function(time_seq)
# Plotting Survival Curve
plot(time_seq, survival_vals, type = "l", col = "blue", lwd = 2,
     xlab = "Time (days)", ylab = "Survival Probability", main = "Survival Curve")

# Plotting Hazard Rate Curve
plot(time_seq, hazard_vals, type = "l", col = "red", lwd = 2,
     xlab = "Time (days)", ylab = "Hazard Rate", main = "Hazard Curve")

# Results
cat("MLE Estimates:\n")
cat("Shape (alpha):", mle_alpha, "\n")
cat("Scale (beta):", mle_beta, "\n")
