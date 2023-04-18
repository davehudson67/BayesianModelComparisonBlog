## exponential
waicEcensus <- function(zL, zU, censored, postSamples) {
  
    ## calculate log likelihood
    log.like <- function(zL, zU, censored, c) {
        ## calculate log-likelihoods interval-censored
        llI <- pexp(zU[censored == 1], c, log.p = TRUE)
        llI <- cbind(llI, pexp(zL[censored == 1], c, log.p = TRUE))
        llI <- apply(llI, 1, function(x) {
            maxx <- max(x)
            x <- exp(x - maxx)
            x <- log(x[1] - x[2])
            maxx + x
        })
        ## calculate log-likelihoods right-censored
        llR <- pexp(zL[censored == 2], c, log.p = TRUE, lower.tail = FALSE)
        l <- c(llI, llR)
    }

    ## calculate likelihoods in parallel
    l <- apply(postSamples, 1, list)
    l <- purrr::map(l, 1)
    l <- mclapply(l,
        function(pars, zL, zU, censored) {
            log.like(zL, zU, censored, pars[1])
        }, zL = zL, zU = zU, censored = censored, mc.cores = 8)
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

## gompertz
waicGcensus <- function(zL, zU, censored, postSamples) {

    ## calculate log likelihood
    log.like <- function(zL, zU, censored, a, b) {
        ## calculate log-likelihoods interval-censored
        llI <- (a / b) * (1 - exp(b * zL[censored == 1]))
        llI <- cbind(llI, (a / b) * (1 - exp(b * zU[censored == 1])))   
        llI <- apply(llI, 1, function(x) {
            maxx <- max(x)
            x <- exp(x - maxx)
            x <- log(x[1] - x[2])
            maxx + x
        })
        ## calculate log-likelihoods right-censored
        llR <- (a / b) * (1 - exp(b * zL[censored == 2]))
        ll <- c(llI, llR)
        ll
    }
  
    ## calculate likelihoods in parallel
    l <- apply(postSamples, 1, list)
    l <- purrr::map(l, 1)
    l <- mclapply(l,
        function(pars, zL, zU, censored) { 
            log.like(zL, zU, censored, pars[1], pars[2])
        }, zL = zL, zU = zU, censored = censored, mc.cores = 8)
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

## gompertz-makeham
waicGMcensus <- function(zL, zU, censored, postSamples) {

    ## calculate log likelihood
    log.like <- function(zL, zU, censored, a, b, c) {
        ## calculate log-likelihoods interval-censored
        llI <- -c * zL[censored == 1] - (a / b) * (exp(b * zL[censored == 1]) - 1)
        llI <- cbind(llI, (-c * zU[censored == 1] - (a / b) * (exp(b * zU[censored == 1]) - 1)))   
        llI <- apply(llI, 1, function(x) {
            maxx <- max(x)
            x <- exp(x - maxx)
            x <- log(x[1] - x[2])
            maxx + x
        })

        ## calculate log-likelihoods right-censored
        llR <- (-c * zL[censored == 2] - (a / b) * (exp(b * zL[censored == 2]) - 1))
        ll <- c(llI, llR)
        ll
    }

    ## calculate likelihoods in parallel
    l <- apply(postSamples, 1, list)
    l <- purrr::map(l, 1)
    l <- mclapply(l,
        function(pars, zL, zU, censored) { 
            log.like(zL, zU, censored, pars[1], pars[2], pars[3])
        }, zL = zL, zU = zU, censored = censored, mc.cores = 8)
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

