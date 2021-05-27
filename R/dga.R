#' Fit decomposable graphical models
#'
#' Fit the decomposable graphical model approach of Madigan and York (1997), as implemented in the dga package of Kristian Lum, James Johndrow, and Patrick Ball (2015).
#'
#' @usage dga(data, priorCounts = NULL, graphs = NULL, Nmissing = NULL,
#'  logPriorGraphs = NULL, logPriorN = NULL, maxMissingMult = 30,
#'  length.out = 1000)
#'
#' @param data object of class "MSEdata" representing list inclusion pattern counts for between 3 and 5 lists.
#' @param priorCounts "Prior counts" for the hyper-Dirichlet prior on decomposable graphs. This can be a positive numeric value for constant counts or an L-dimensional array (where L is the number of lists) representing the table of prior counts (see examples below).
#' @param graphs List of decomposable graphs under consideration. Must be a sub-list of "graphs3", "graphs4" or "graphs5". Default is "graphsL" where L is the number of lists (between 3 and 5).
#' @param Nmissing integer vector for the number of plausible unobserved individuals. Population size posterior probability is computed for the values corresponding to these numbers of missing individuals. Defaults to a vector of length 1000 ranging between 1 and 30 times the number of observed individuals.
#' @param logPriorGraphs Numeric vector of log prior probabilities over the set of decomposable graphs in "graphs". Defaults to a constant prior.
#' @param logPriorN Numeric vector of log prior probabilities of the number of missing individuals in "Nmissing". Defaults to the population size prior $p(N) \propto 1/N$. Nmissing should be specified explicitely if logPriorN is provided.
#' @param maxMissingMult If Nmissing is NULL, then maxMissingMult determines the default range of Nmissing: between 1 and maxMissingMult times the number of observed individuals.
#' @param length.out If Nmissing is NULL, then length.out determines the maximum default length of Nmissing.
#'
#' @value Object of class "dga" containing the population sizes "N", the "weight" matrix of posterior probabilities for each population size and graphical model, and the "args" list of arguments passed to this function. Use the \code{estimates()} function to obtain point and interval estimates.
#'
#' @examples
#' dga_fit <- dga(UK)
#' estimates(dga_fit)
#' plot(dga_fit)
#'
#' @seealso estimates MSEdata
#'
#' @import assert dga
#' @export
dga <- function(data,
                priorCounts = NULL,
                graphs = NULL,
                Nmissing = NULL,
                logPriorGraphs = NULL,
                logPriorN = NULL,
                maxMissingMult = 30,
                length.out = 1000) {
  assert(is.MSEdata(data))
  args <- as.list(environment())

  if (is.null(priorCounts)) {
    priorCounts <- 2^(-nlists(data))
  }

  if (is.null(graphs)) {
    if (nlists(data) == 3) {
      load(system.file("data/graphs3.rda", package="dga"))
      graphs <- graphs3
    } else if (nlists(data) == 4) {
      load(system.file("data/graphs4.rda", package="dga"))
      graphs <- graphs4
    } else if (nlists(data) == 5) {
      load(system.file("data/graphs5.rda", package="dga"))
      graphs <- graphs5
    }
  }

  if (is.null(Nmissing)) {
    Nmissing <- seq(1, maxMissingMult * sum(data$count), length.out = length.out)
  }
  Nmissing <- unique(round(Nmissing))

  if (is.null(logPriorN)) {
    logPriorN <- -log(sum(data$count) + Nmissing)
  }

  Y <- MSEdata_to_array(data)
  # Fit dga
  weights <- bma.cr(Y, Nmissing, priorCounts, graphs,
    logprior = logPriorN,
    log.prior.model.weights = logPriorGraphs
  )

  # Construct dga object
  dga_fit <- list(
    weights = weights,
    N = Nmissing + sum(data$count),
    args = args
  )
  structure(dga_fit, class = c("dga", "MSEfit"))
}

is.dga <- function(x) {
  inherits(x, "dga")
}
