library(loo)

## load data and samples
load("Sim_Data.RData")
mcmc.out <- readRDS("samples/mcmc.out_h1.rds")
mcmc_samples <- as.matrix(mcmc.out$samples)
mcmc_samples <- mcmc_samples[1:1000,]

## define log likelihood function
log.like <- function(x, y, beta0, beta1, sigma) {
  mu <- beta0 + beta1 * x
  ll <- dnorm(y, mu, sigma, log = TRUE)
  return(ll)
}

## calculate likelihoods in parallel from MCMC samples and data
l <- apply(mcmc_samples, 1, list)
l <- purrr::map(l, 1)
l <- mclapply(l,
              function(pars, x, y) { 
                log.like(x, y, pars[1], pars[2], pars[3])
              }, x = x, y = y, mc.cores = 8)

## create matrix of likelihoods (columns = individuals, rows = samples)
l <- reduce(l, base::rbind)

## create vector to identify which chain each set of samples came from
chain_id <- rep(1:2, each = 500)

## calculate the MCMC effective samples size divided by the total sample size
reff <- relative_eff(l, chain_id)

## compute PSIS-LOO CV
loo(l, r_eff = reff)

## calculate WAIC using loo package
WAIC_loo <- loo::waic(l)
WAIC_loo_h1 <- WAIC_loo$estimates[3]
