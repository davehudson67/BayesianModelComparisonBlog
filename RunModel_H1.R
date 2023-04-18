## set up Nimble code
code <- nimbleCode({
  
  ## likelihood
  for (i in 1:n){
    y[i] ~ dnorm(mu[i], sd = sigma)
    mu[i] <- beta0 + beta1 * x[i]
  }
  
  ## priors
  beta0 ~ dnorm(0, 0.001)
  beta1 ~ dnorm(0, 0.001)
  sigma ~ dunif(0, 100)

})

consts <- list(n = length(x))
data <- list(y = y, x = x)
inits <- list(beta0 = rnorm(1),
              beta1 = rnorm(1),
              sigma = runif(1, 0, 100))

model <- nimbleModel(code = code, constants = consts, data = data, inits = inits)

cmodel <- compileNimble(model)

config <- configureMCMC(cmodel, monitors = c("beta0", "beta1", "sigma"), thin = 1, enableWAIC = TRUE)
config$removeSamplers(c("beta0", "beta1"))
config$addSampler(target = c("beta0", "beta1"), type = "AF_slice")


built <- buildMCMC(config)

cBuilt <- compileNimble(built)

mcmc.out <- runMCMC(cBuilt,
                    niter = 100000,
                    summary = TRUE,
                    WAIC = TRUE,
                    samplesAsCodaMCMC = TRUE,
                    progressBar = TRUE,
                    nchains = 2,
                    nburnin = 1000)

saveRDS(mcmc.out, "samples/mcmc.out_h1.rds")

mcmc.out$WAIC
mcmc.out$summary
plot(mcmc.out$samples)
WAIC_nimble_h1 <- mcmc.out$WAIC$WAIC
