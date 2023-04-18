## load data
load("Sim_Data.RData")
mcmc.out <- readRDS("samples/mcmc.out_h0.rds")
mcmc_samples <- as.matrix(mcmc.out$samples)
mcmc_samples <- mcmc_samples[1:1000,]


waic <- function(y, postSamples) {
  
  log.like <- function(y, beta0, sigma) {
    mu <- beta0
    ll <- dnorm(y, mu, sigma, log = TRUE)
    return(ll)
  }
  
  ## calculate likelihoods in parallel
  l <- apply(postSamples, 1, list)
  l <- purrr::map(l, 1)
  l <- mclapply(l,
                function(pars, y) { 
                  log.like(y, pars[1], pars[2])
                }, y = y, mc.cores = 8)
  ## create matrix of likelihoods (columns = individuals, rows = samples)
  l <- reduce(l, base::rbind)
  
  ## now deal with first half of the equation - mean of the likelihood for each individual observation across the samples
  lppd <- apply(l, 2, log_sum_exp_marg, mn = TRUE)
  lppd <- sum(lppd)  
  
  ## now deal with the second half of the equation - mean sum of the log likelihood for each individual across the samples
  l2 <- colMeans(l)
  l2 <- sum(l2)
  
  pwaic1 <- 2 * (lppd - l2)
  
  waic <- -2 * (lppd - pwaic1)
  
  return(waic)
}

WAIC_naive_h0 <- waic(y, mcmc_samples)
