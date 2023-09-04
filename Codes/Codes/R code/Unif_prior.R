library(mcmc)
library(gtools)
#Load the whole sample
sample<-read.csv('/Users/shambhusharansinha/Downloads/data2 (1).csv')
# Get the size of the array
array_size <- length(sample)

# Print the size of the array##
cat("Size of the array:", array_size, "\n")
data<-sample[1:10000,1]
log_likelihood <- function(a, b, c, data1) {
  m <- length(data1)
  l <- (m * b) / a + m * log(a) - m * c * log(a)+m*c*log(b)-m*log(gamma(c))
  for (i in 1:m) {
    if(exp(a * data1[i]) - 1<=0)
    {
      return -100000
    }
    l <- l + a * data1[i] + (c - 1) * (log(exp(a * data1[i]) - 1)) - (b / a) * exp(a * data1[i])
  }
  return(l)
}

# Define log-prior function
log_prior <- function(a, b, c) {
  if(a<0||a>10||b<0||b>10||c<0||c>10)
  {
    log_prior<- -1000000
  }
  log_prior<- -3 #uniform on (0,10))
}

# Define log-posterior function
log_posterior <- function(a, b, c, data) {
  log_likelihood(a, b, c, data) + log_prior(a, b, c)
}


# Set up MCMC parameters
nsteps <- 100000
step_size <- c(0.01, 0.01, 0.01)
init_params <- c(1.1,1.8,0.7)
samples <- matrix(0, nrow = nsteps, ncol = 3)
current_params <- init_params

# Run MCMC sampling
for (step in 1:nsteps) {
  # Propose new parameters
  proposed_params <- current_params + step_size * rnorm(3)
  # Check if b is not equal to 0
  
  if(proposed_params[2]>0 && proposed_params[1]>0)
  { # Compute log-ratio of posterior probabilities
    log_ratio <- log_posterior(proposed_params[1], proposed_params[2], proposed_params[3], data) -
      log_posterior(current_params[1], current_params[2], current_params[3], data)
    
    # Accept or reject the proposed parameters
    if (log(runif(1)) < log_ratio) {
      current_params <- proposed_params
      
    }
  }
  
  # Store samples
  if(step>1000)
  {samples[step, ] <- current_params}
}

# Plot the trace of samples
matplot(samples, type = 'l', xlab = 'Step', ylab = 'Parameter Value', main = 'MCMC Trace')
legend('topright', c('a', 'b', 'c'), lty = 1)

# Compute posterior mean and posterior variance of samples
posterior_mean <- colMeans(samples)
posterior_variance <- colMeans((samples - posterior_mean)^2)
posterior_mean
posterior_variance

# Compute credible intervals (2.5th and 97.5th percentiles)
credible_interval <- function(samples, alpha) {
  n <- nrow(samples)
  lower <- floor(n * alpha / 2)
  upper <- ceiling(n * (1 - alpha / 2))
  cred_int <- apply(samples, 2, function(x) quantile(x, probs = c(alpha/2, 1 - alpha/2)))
  return(cred_int)
}

# Usage
alpha <- 0.05 # alpha level for credible intervals (e.g., 95% CI corresponds to alpha = 0.05)
cred_int <- credible_interval(samples, alpha)
cat("Credible Intervals (", 100 * (1 - alpha), "%):", "\n")
cat("Parameter a: [", cred_int[1, 1], ", ", cred_int[2, 1], "]\n", sep = "")
cat("Parameter b: [", cred_int[1, 2], ", ", cred_int[2, 2], "]\n", sep = "")
cat("Parameter c: [", cred_int[1, 3], ", ", cred_int[2, 3], "]\n", sep = "")

alpha <- 0.01 # alpha level for credible intervals (e.g., 95% CI corresponds to alpha = 0.05)
cred_int <- credible_interval(samples, alpha)
cat("Credible Intervals (", 100 * (1 - alpha), "%):", "\n")
cat("Parameter a: [", cred_int[1, 1], ", ", cred_int[2, 1], "]\n", sep = "")
cat("Parameter b: [", cred_int[1, 2], ", ", cred_int[2, 2], "]\n", sep = "")
cat("Parameter c: [", cred_int[1, 3], ", ", cred_int[2, 3], "]\n", sep = "")


