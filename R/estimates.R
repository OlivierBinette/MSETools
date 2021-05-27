#'
#' @export
estimates <- function(fit, ...) {
  UseMethod("estimates", fit)
}

#' @export
estimates.list <- function(fit, ...) {
  lapply(fit, function(x) estimates(x, ...))
}

#' @export
estimates.independence <- function(fit, level = 0.95) {
  assert(is.independence(fit))

  nobs = sum(fit$args$data$count)

  est <- exp(coef(fit$fit$fit)[["(Intercept)"]]) + nobs
  CI <- suppressMessages(exp(confint(fit$fit$fit, "(Intercept)"))) + nobs

  return(c(
    N.hat = est,
    N.lwr = CI[[1]],
    N.upr = CI[[2]]
  ))
}

#' @export
estimates.dga <- function(fit, level = 0.95) {
  assert(is.dga(fit))

  weights <- colSums(fit$weights)
  N <- fit$N

  cs <- cumsum(weights)
  median <- N[which(cs > 0.5)[1]]

  levels <- c((1 - level) / 2, (1 + level) / 2)
  CI <- sapply(levels, function(level) N[which(cs > level)[1]])

  return(c(
    N.hat = median,
    N.lwr = CI[[1]],
    N.upr = CI[[2]]
  ))
}

#' @import parallel
#' @export
estimates.lcmcr <- function(fit, level = 0.95, burnin = 0, nSamples = 100, thinning = 1000, mc.cores=detectCores()) {
  assert(is.lcmcr(fit))

  # Always re-run.
  args = fit$args
  args$lazy=FALSE
  fit <- do.call(lcmcr, args)

  samples <- mclapply(fit$lcmcr_chains, function(lcmcr) {
    lcmcr$Update(burnin, output = FALSE)
    lcmcr$Change_Trace_Length(nSamples)
    lcmcr$Change_SubSamp(thinning)
    lcmcr$Activate_Tracing()

    # Sample
    lcmcr$Update(nSamples * thinning, output = FALSE)

    return(lcmcr$Get_Trace("n0") + lcmcr$n)
  }, mc.cores = mc.cores, mc.preschedule = FALSE)

  samples <- unlist(samples)
  CI <- quantile(samples, c((1 - level) / 2, 0.5, (1 + level) / 2))

  ests <- c(CI[2], CI[1], CI[3])
  names(ests) <- c("N.hat", "N.lwr", "N.upr")

  return(ests)
}

#' @import parallel
#' @export
estimates.sparsemse <- function(fit, level = 0.95, nboot = 2000, iseed = 1234, mc.cores = detectCores()) {
  assert(is.sparsemse(fit))

  args <- fit$args
  data <- args$data
  pthresh <- args$pthresh
  alpha <- c((1 - level) / 2, (1 + level) / 2)

  # Parallelized adaptation of SparseMSE::estimatepopulation
  RNGkind("L'Ecuyer-CMRG")
  set.seed(iseed)
  n1 <- dim(data)[1]
  n2 <- dim(data)[2]
  countsobserved <- data[, n2]
  nobs <- sum(countsobserved)
  populationestimatefromdata <- fit$fit

  popest <- populationestimatefromdata$estimate
  MSEfit <- populationestimatefromdata$MSEfit

  bootreps <- parallel::mcmapply(function(i) {
    counts <- rmultinom(1, nobs, countsobserved)
    zdatboot <- cbind(data[, -n2], counts)
    SparseMSE::estimatepopulation.0(zdatboot,
                                    quantiles = NULL,
                                    pthresh = pthresh
    )$estimate
  }, 1:nboot, mc.preschedule = TRUE, mc.cores = mc.cores)

  # Reproduced from SparseMSE::estimatepopulation
  jackest <- rep(0, n1)
  for (j in (1:n1)) {
    nj <- data[j, n2]
    if (nj > 0) {
      zd1 <- data
      zd1[j, n2] <- nj - 1
      jackest[j] <- SparseMSE::estimatepopulation.0(zd1,
                                                    quantiles = NULL,
                                                    pthresh = pthresh
      )$estimate
    }
  }
  jr <- sum(countsobserved * jackest) / sum(countsobserved) -
    jackest
  ahat <- sum(countsobserved * jr^3) / (6 * (sum(countsobserved * jr^2))^{
    3 / 2
  })
  confquantiles <- SparseMSE::bcaconfvalues(bootreps, popest, ahat, alpha)

  ests <- c(popest, confquantiles[[1]], confquantiles[[2]])
  names(ests) <- c("N.hat", "N.lwr", "N.upr")

  return(ests)
}
