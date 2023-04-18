##Functions used in Model Comparisons:

## calculates log of the marginal likelihood using log-sum-exp trick
log_sum_exp_marg <- function(x, ind = NULL, mn = TRUE) {
  
  if(!is.null(ind)) {
    x <- x[ind]
  }
  
  ## extract length
  n <- length(x)
  
  ## extract non-finite values
  ## (it's OK to remove these here because we
  ## assume they correspond to a likelihood of zero
  ## and are thus felt through 1 / n above
  ## in the marginal likelihood calculation)
  x <- x[is.finite(x)]
  if(length(x) == 0) {
    return(-Inf)
  }
  
  ## extract maximum of logged values
  mX <- max(x)
  
  ## return answer
  out <- mX + log(sum(exp(x - mX)))
  if(mn) {
    out <- out - log(n)
  }
  out
}

## create and plot bootstrap samples of importance weights with 95% CIs
BootsPlot <- function(impWeights, r, trace = FALSE, nsteps = 10) {
  #browser()
  if(!trace) {
      boot.iw <- boot(data = impWeights, statistic = log_sum_exp_marg, R = r)
      ci <- boot.ci(boot.iw, type = "basic")
      hist(boot.iw$t)
      abline(v=ci$basic[4:5], col="red")
      br <- list(logmarg=boot.iw$t0, LCI=ci$basic[4],UCI=ci$basic[5])
  } else {
    
    ntot <- length(impWeights)
    if(ntot < 10) {
      stop("Number of importance samples < 10")
    }
    nsamps <- seq(10, ntot, length.out = nsteps)
    nsamps <- unique(round(nsamps))
    out <- NULL
    for (i in 1:length(nsamps)) {
      br <- impWeights[sample(length(impWeights), nsamps[i], replace = FALSE)]
      br <- boot(data = br, statistic = log_sum_exp_marg, R = r)
      br <- boot.ci(br, type = "basic")
      br <- data.frame(logmarg = br$t0, LCI = br$basic[4], UCI = br$basic[5])
      out[[i]] <- br
    }
    names(out) <- nsamps
    out <- bind_rows(out, .id = "nsamp") %>%
      mutate(nsamp = as.numeric(nsamp))
    
    print(ggplot(out, aes(x = nsamp)) +
      geom_point(aes(y = logmarg)) +
      geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
      geom_hline(yintercept = out$logmarg[out$nsamp == max(out$nsamp)], linetype = "dashed"))
  }
  br
}

## create marginal likelihood plots
MargLike.plot <- function(boots, thresh = log(20)) {
  m <- purrr::map(boots, as.data.frame)
  br <- bind_rows(m, .id = "model")
  br$model <- as.factor(br$model)
  
  p <- ggplot(br, aes(x=model, y=logmarg)) +
    geom_point(aes(size=5)) +
    geom_errorbar(aes(ymin=LCI, ymax=UCI)) +
    labs(x = "Model", y = "Log marginal likelihood") +
#    theme_bw() +
#    scale_x_discrete(labels = c("No Escort differences", "Escort difference a1", "Escort difference a2", "Escort difference b1", "Escort difference b2", "Escort difference c1")) +
    theme(legend.position = "none", axis.text = element_text(size = 12), axis.title = element_text(size = 14))
  
  if(!is.na(thresh)) {
      ## which models within thresh of best
      logmarg <- map_dbl(mods, "logmarg")
      bestind <- which(logmarg == max(logmarg))
      logmargLCI <- mods[[bestind]]$LCI
      thresh <- logmargLCI - thresh
      p <- p + geom_hline(yintercept = thresh, linetype = "dashed")
  }
  p
} 
  




