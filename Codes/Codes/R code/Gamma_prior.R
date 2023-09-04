library(mcmc)
library(boot)
library(gtools)

# Load the whole sample
sample <- read.csv('/Users/shambhusharansinha/Downloads/data1 (1).csv')
# Sample 1000 random rows from the data frame
random_indices <- sample(nrow(sample), 10000, replace = FALSE)
data <- data[random_indices, ]
data <- as.data.frame(data)
log_likelihood <- function(a, b, c, data1) {
  m <- length(data1)
  l <- (m * b) / a + m * log(a) - m * c * log(a) + m * c * log(b) - m * log(gamma(c))
  for (i in 1:m) {
    if (exp(a * data1[i]) - 1 <= 0) 
      {
      return -100000
    }
    l <- l + a * data1[i] + (c - 1) * (log(exp(a * data1[i]) - 1)) - (b / a) * exp(a * data1[i])
  }
  return(l)
}

# Define log-prior function
log_prior <- function(a, b, c) {
  sum(dgamma(a, shape = 1, rate = 1, log = TRUE)) +
    sum(dgamma(b, shape = 1, rate = 1, log = TRUE)) +
    sum(dgamma(c, shape = 1, rate = 1, log = TRUE))
}

# Define log-posterior function
log_posterior <- function(a, b, c, data) {
  log_likelihood(a, b, c, data) + log_prior(a, b, c)
}

# Set up MCMC parameters
nsteps <- 10
step_size <- c(0.01, 0.01, 0.01)
init_params <- c(1.4, 2.2, 0.8)
samples <- matrix(0, nrow = nsteps, ncol = 3)
current_params <- init_params

# Run MCMC sampling
for (step in 1:nsteps) {
  # Propose new parameters
  proposed_params <- current_params + step_size * rnorm(3)
  # Check if b is not equal to 0
  if (proposed_params[2] > 0 && proposed_params[1] > 0) {
    # Compute log-ratio of posterior probabilities
    log_ratio <- log_posterior(proposed_params[1], proposed_params[2], proposed_params[3], data) -
      log_posterior(current_params[1], current_params[2], current_params[3], data)
    
    # Accept or reject the proposed parameters
    if (log(runif(1)) < log_ratio) {
      current_params <- proposed_params
    }
  }
  
  # Store samples
  if (step > 1) {
    samples[step, ] <- current_params
  }
}
boot_func <- function(data, indices) {
  boot_data <- data[indices, ]  # Subset data using indices
  boot_samples <- matrix(0, nrow = nsteps, ncol = 3)
  boot_current_params <- init_params
  for (step in 1:nsteps) {
    # Propose new parameters
    proposed_params <- boot_current_params + step_size * rnorm(3)
    # Check if b is not equal to 0
    if (proposed_params[2] > 0 && proposed_params[1] > 0) {
      # Compute log-ratio of posterior probabilities
      log_ratio <- log_posterior(proposed_params[1], proposed_params[2], proposed_params[3], boot_data) -
        log_posterior(boot_current_params[1], boot_current_params[2], boot_current_params[3], boot_data)
      
      # Accept or reject the proposed parameters
      if (log(runif(1)) < log_ratio) {
        boot_current_params <- proposed_params
      }
    }
    
    # Store samples
    if (step > 1) {
      # Randomly sample 1000 indices without replacement
      sample_indices <- sample(1:step, 1000, replace = FALSE)
      boot_samples[sample_indices, ] <- boot_current_params
    }
  }
  
  return(boot_samples)
}

# Set seed for reproducibility
set.seed(123)

# Run bootstrapping
boot_results <- boot(data, statistic = boot_func, R = 1)
# Print estimated values and confidence intervals
cat("Estimated values:\n")
print(colMeans(boot_results$t))

cat("\nVariance:\n")
print(var(boot_results$t, na.rm = TRUE))

cat("\n95% Confidence intervals:\n")
boot_ci_95 <- boot.ci(boot_results, type = "bca", index = 1:3, conf = 0.95)
print(boot_ci_95$percent[4:6, ])

cat("\n97% Confidence intervals:\n")
boot_ci_97 <- boot.ci(boot_results, type = "bca", index = 1:3, conf = 0.97)
print(boot_ci_97$percent[4:6, ])

cat("\n99% Confidence intervals:\n")
boot_ci_99 <- boot.ci(boot_results, type = "bca", index = 1:3, conf = 0.99)
print(boot_ci_99$percent[4:6, ])
