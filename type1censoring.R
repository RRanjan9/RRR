#Initial information
n <- 500              # Total number of units
true_lambda <- 0.5    # The actual rate parameter for the exponential distribution
censor_time <- 5      # Time at which we stop observing (Type 1 censoring)

# Generating failure times using the exponential distribution
failure_times <- rexp(n, rate = true_lambda)  

# Applying censoring
observed_times <- pmin(failure_times, censor_time)  # Take the minimum of failure time or censor time
censored <- failure_times > censor_time             # Check which times are censored (TRUE/FALSE)

# Counting failures and censored cases
m <- sum(!censored)            # Number of units that failed (uncensored)
n_censored <- sum(censored)    # Number of units that are censored (still functioning at censor time)

# Step 5: Calculate the estimate of lambda
m <- sum(observed_times[!censored])  # Total failure time for units that failed
lambda_hat <- m / (m + n_censored * censor_time)  # Formula to calculate estimated lambda

# Step 6: Show the results
cat("True lambda:", true_lambda, "\n")
cat("Estimated lambda:", lambda_hat, "\n")
cat("Number of failures:", r, "\n")
cat("Number of censored observations:", n_censored, "\n")

# Load necessary libraries
library(ggplot2)

# Create a data frame with the observed times and censoring status
data <- data.frame(
  observed_times = observed_times,
  status = ifelse(censored, "Censored", "Failed")
)

# Generate a histogram to visualize the distribution of failure and censored times
ggplot(data, aes(x = observed_times, fill = status)) +
  geom_histogram(binwidth = 0.5, alpha = 0.7, position = "identity", color = "black") +
  scale_fill_manual(values = c("Censored" = "red", "Failed" = "blue")) +
  labs(title = "Distribution of Failure Times with Type 1 Censoring",
       x = "Observed Times", y = "Count", fill = "Status") +
  geom_vline(xintercept = censor_time, linetype = "dashed", color = "black") +
  theme_minimal()

