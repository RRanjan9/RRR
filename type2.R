# Set parameters for the exponential distribution
set.seed(123)
n <- 100           # Total number of units
true_lambda <- 0.3 # True rate parameter (1/lambda is the mean)
censor_time <- 5   # Type 1 censoring time

# Generate n samples from an exponential distribution
failure_times <- rexp(n, rate = true_lambda)

# Apply Type 1 censoring
censored <- failure_times > censor_time
observed_times <- pmin(failure_times, censor_time) # Take min of failure time and censor time

# Count the number of failures and censored units
r <- sum(!censored) # Number of failures (uncensored)
n_censored <- sum(censored) # Number of censored units

# Calculate the MLE of lambda under Type 1 censoring
T <- sum(observed_times[!censored]) # Sum of failure times (uncensored)
lambda_hat <- r / (T + n_censored * censor_time) # MLE of lambda

# Output the results
cat("True lambda:", true_lambda, "\n")
cat("Estimated lambda:", lambda_hat, "\n")
cat("Number of failures:", r, "\n")
cat("Number of censored observations:", n_censored, "\n")
