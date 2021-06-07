
#' Get MCMC traces for LCMCR
#'
#' Return MCMC traces for the "n0", "prob_zero" and "k_star" parameters of LCMCR objects. For each of these parameters, a matrix is returned where each row represents an iteration and each column a chain.
#'
#' @usage MCMCtraces(fit, burnin = 0, nSamples = 100, thinning = 1000, mc.cores=detectCores())
#'
#' @import purrr dplyr parallel
#' @export
MCMCtrace <- function(fit, burnin = 0, nSamples = 100, thinning = 1000, mc.cores=detectCores()) {
  assert(is.lcmcr(fit))

  # Always re-run.
  fit <- lcmcr(fit$args$data, K = fit$args$K, seeds = fit$args$seeds, lazy=FALSE)
  traces <- mclapply(fit$lcmcr_chains, function(lcmcr) {
    lcmcr$Update(burnin, output = FALSE)
    lcmcr$Change_Trace_Length(nSamples)
    lcmcr$Change_SubSamp(thinning)
    lcmcr$Activate_Tracing()

    # Sample
    lcmcr$Update(nSamples * thinning, output = FALSE)

    trace = list(n0 = lcmcr$Get_Trace("n0"),
                 prob_zero = lcmcr$Get_Trace("prob_zero"),
                 k_star = lcmcr$Get_Trace("k_star"))
    return(trace)
  }, mc.cores = mc.cores, mc.preschedule = FALSE)

  n0 = sapply(map(traces, "n0"), function(x) x)
  prob_zero = sapply(map(traces, "prob_zero"), function(x) x)
  k_star = sapply(map(traces, "k_star"), function(x) x)

  return(list(n0 = n0,
              prob_zero = prob_zero,
              k_star = k_star))
}

#' @importFrom rstan Rhat ess_bulk
#' @export
diagnostics <- function(traces) {

  functs = list(Rhat = rstan::Rhat,
                n_ess = rstan::ess_bulk,
                n_iter = function(x) length(x))

  lapply(functs, function(FUN) {
    sapply(traces, function(trace) FUN(trace))
  }) %>%
    data.frame()
}
