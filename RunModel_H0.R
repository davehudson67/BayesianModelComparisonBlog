## set up Nimble code
code <- nimbleCode({
  
  ## likelihood
  for (i in 1:n){
    y[i] ~ dnorm(mu[i], sd = sigma)
    mu[i] <- beta0
  }
  
  ## priors
  beta0 ~ dnorm(0, sd = 1)
  sigma ~ dunif(0, 100)
  
})

consts <- list(n = length(x))
data <- list(y = y)
inits <- list(beta0 = rnorm(1),
              sigma = runif(1, 0, 100))

model <- nimbleModel(code = code, constants = consts, data = data, inits = inits)

cmodel <- compileNimble(model)

config <- configureMCMC(cmodel, monitors = c("beta0", "sigma"), thin = 1, enableWAIC = TRUE)
config$removeSamplers(c("beta0"))
config$addSampler(target = c("beta0"), type = "slice")

built <- buildMCMC(config)

cBuilt <- compileNimble(built)

mcmc.out <- runMCMC(cBuilt,
                    niter = 10000,
                    summary = TRUE,
                    WAIC = TRUE,
                    samplesAsCodaMCMC = TRUE,
                    progressBar = TRUE,
                    nchains = 2,
                    nburnin = 1000)

saveRDS(mcmc.out, "samples/mcmc.out_h0.rds")

mcmc.out$WAIC
mcmc.out$summary
plot(mcmc.out$samples)
WAIC_nimble_h0 <- mcmc.out$WAIC$WAIC
