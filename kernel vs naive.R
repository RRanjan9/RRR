# Dataset
data <- c(5.65746599, 5.38283914, 2.79892121, 2.85423660, 2.95252721, 5.42626667,
          7.66239113, -0.18001073, 0.65083500, 2.40276530, -0.09929884, 6.32619215,
          5.03650752, 2.07470777, 1.78019174, 6.12891558, 4.05352439, 2.02686971,
          3.50834853, -2.76449768, 4.98428763, 3.01292677, 2.82448038, 3.98110437,
          5.09371862, 5.97961648, 4.56968496, -0.48814532, 5.08736697, 2.41757609)

# Kernel Density Estimate
kde <- density(data)
kde

# Naïve Density Estimator
naive_density <- function(x, data, h) {
  n <- length(data)
  return(sum((abs(x - data) <= h) / (2 * h)) / n)
}

# Naïve Density Estimator for a range of x
h <- 1.0 # Bandwidth for Naïve Density Estimation
x_range <- seq(min(data) - 1, max(data) + 1, length.out = 100)
naive_estimates <- sapply(x_range, naive_density, data = data, h = h)

# Plotting both
plot(kde, main = "Kernel Density Estimate vs Naïve Density Estimator", 
     col = "blue", lwd = 2, ylim = c(0, max(kde$y, naive_estimates)), xlab = "X", ylab = "Density")
lines(x_range, naive_estimates, col = "red", lwd = 2)
legend("topright", legend = c("KDE", "Naïve Density"), col = c("blue", "red"), lty = 1, lwd = 2)
