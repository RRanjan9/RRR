shape_param <- 3     # Shape 
scale_param <- 10    # Scale
n_samples <- 100     # Number of observations

# Weibull-distributed data
data <- scale_param * (-log(runif(n_samples)))^(1 / shape_param)

# Negative log-likelihood function
neg_log_likelihood <- function(alpha, beta) {
  if (alpha <= 0 || beta <= 0) return(Inf) # Ensure parameters are positive
  n <- length(data)
  log_likelihood <- 
    n * log(alpha) - n * log(beta) +
    (alpha - 1) * sum(log(data)) -
    sum((data / beta)^alpha)
  return(-log_likelihood)  # Negative log-likelihood
}

# Generate 2D grid for shape (alpha) and scale (beta)
alpha_vals <- seq(2, 4, length.out = 50)  # Shape parameter grid
beta_vals <- seq(8, 12, length.out = 50)  # Scale parameter grid
likelihood_matrix <- matrix(0, nrow = length(alpha_vals), ncol = length(beta_vals))

# Computing log-likelihood for each combination
for (i in 1:length(alpha_vals)) {
  for (j in 1:length(beta_vals)) {
    likelihood_matrix[i, j] <- -neg_log_likelihood(alpha_vals[i], beta_vals[j])
  }}
# 2D likelihood plot
filled.contour(
  x = alpha_vals, y = beta_vals, z = likelihood_matrix,
  xlab = "Shape (alpha)", ylab = "Scale (beta)",
  main = "2D Likelihood Plot of Weibull Model"
)
# MLE Estimation
result <- optim(par = c(2, 5), fn = function(params) neg_log_likelihood(params[1], params[2]), 
                method = "L-BFGS-B", lower = c(1e-5, 1e-5))
mle_shape <- result$par[1]
mle_scale <- result$par[2]

# Mean failure time
mean_failure_time <- mle_scale * gamma(1 + 1 / mle_shape)  # Using Weibull formula
sample_mean <- mean(data)  # Sample mean
#results
cat("MLE Shape (alpha):", mle_shape, "\n")
cat("MLE Scale (beta):", mle_scale, "\n")
cat("MLE Mean Failure Time:", mean_failure_time, "\n")
cat("Sample Mean Failure Time:", sample_mean, "\n")
