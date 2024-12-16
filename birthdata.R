data("birthwt", package = "MASS")
birthwt <- birthwt
birthwt$low <- as.numeric(birthwt$low)  


# Split data into training and test sets
train_indices <- sample(1:nrow(birthwt), size = 0.7 * nrow(birthwt))  # 70% training
train_data <- birthwt[train_indices, ]
test_data <- birthwt[-train_indices, ]

# Logistic function
logistic <- function(x) {
  1 / (1 + exp(-x))
}

# Fit logistic regression using Newton-Raphson
fit_logistic <- function(X, y, max_iter = 100, tol = 1e-6) {
  n <- nrow(X)
  p <- ncol(X)
  
  # Initialize beta coefficients
  beta <- matrix(0, nrow = p, ncol = 1)
  epsilon <- 1e-8  
  
  for (i in 1:max_iter) {
    eta <- X %*% beta  # Linear predictor
    mu <- logistic(eta)  # Predicted probabilities
    mu <- pmin(pmax(mu, epsilon), 1 - epsilon)
    
    # Compute gradient and Hessian
    W <- diag(as.vector(mu * (1 - mu)), n, n)  
    z <- eta + (y - mu) / (mu * (1 - mu))
    
    Hessian <- t(X) %*% W %*% X
    gradient <- t(X) %*% (y - mu)
    beta_new <- beta + solve(Hessian + diag(epsilon, p)) %*% gradient  # Regularized update
    
    if (any(is.na(beta_new))) stop("Divergence detected in Newton-Raphson method.")
    if (max(abs(beta_new - beta)) < tol) break
    
    beta <- beta_new
  }
  return(beta)
}

# Addding intercept to training data
X_train <- cbind(1, as.matrix(train_data[, -1]))  # Intercept and regressors
y_train <- train_data$low

# Fitting logistic regression model
beta <- fit_logistic(X_train, y_train)

# Logistic model equation
logistic_model_eq <- paste("log(odds) =", paste(round(beta, 4), c("(Intercept)", colnames(train_data)[-1]), collapse = " + "))
cat("Logistic Model Equation:\n", logistic_model_eq, "\n\n")

# Predict on test data
X_test <- cbind(1, as.matrix(test_data[, -1]))
y_test <- test_data$low
log_odds <- X_test %*% beta
probs <- logistic(log_odds)
predictions <- ifelse(probs > 0.5, 1, 0)

# Calculating confusion matrix
confusion_matrix <- table(Predicted = predictions, Actual = y_test)
cat("Confusion Matrix:\n")
print(confusion_matrix)

# Extracting metrics from confusion matrix
true_positive <- confusion_matrix["1", "1"]
true_negative <- confusion_matrix["0", "0"]
false_positive <- confusion_matrix["1", "0"]
false_negative <- confusion_matrix["0", "1"]

sensitivity <- true_positive / (true_positive + false_negative)
specificity <- true_negative / (true_negative + false_positive)
positive_predictive_value <- true_positive / (true_positive + false_positive)
negative_predictive_value <- true_negative / (true_negative + false_negative)
# Printing evolution metrics
cat("\nModel Performance Metrics:\n")
cat("Sensitivity (True Positive Rate):", round(sensitivity, 4), "\n")
cat("Specificity (True Negative Rate):", round(specificity, 4), "\n")
cat("Positive Predictive Value (PPV):", round(positive_predictive_value, 4), "\n")
cat("Negative Predictive Value (NPV):", round(negative_predictive_value, 4), "\n")
