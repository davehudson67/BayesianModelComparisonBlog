# load MCMC samples for the model parameters
library(coda)

mcmc_samples <- as.matrix(mcmc.out$samples)

# define the likelihood function
likelihood <- function(beta, sigma, y) {
  mu <- beta
  dnorm(y, mu, sigma, log = TRUE)
}

# calculate the deviance of the model for each MCMC sample
deviances <- sapply(1:nrow(mcmc_samples), function(i) {
  #browser()
  beta <- mcmc_samples[i, 1]
  sigma <- mcmc_samples[i, 2]
  -2 * sum(likelihood(beta, sigma, y))
})

# calculate the posterior mean deviance
post_mean_deviance <- mean(deviances)

# calculate the DICp for each parameter
deviance_int <- post_mean_deviance - deviances
p_D <- sum(deviance_int)

# calculate the DIC
DIC <- post_mean_deviance + p_D

