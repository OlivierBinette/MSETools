#' @title Multiple Systems Estimation Using Decomposable Graphical Models
#'
#' @description Compute population size posterior distributions for decomposable graphical models.
#'
#' @usage bma.cr(Y, Nmissing, delta, graphs,
#'               logprior = NULL, log.prior.model.weights = NULL, alpha=1)
#'
#' @param Y \code{p}-dimensional array (\code{2^p} elements) of list intersection counts.
#' @param Nmissing Vector of all possible values for the number of individuals that appear on no list.
#' @param delta Hyper-parameter for the hyper-Dirichlet prior distribution on list intersection probabilities. Can be a positive number of an array of "prior counts" of the same dimension of `Y`.
#' @param graphs Pre-computed list of all decomposable graphical models for \code{p} lists. These should be loaded using data(graphs3), data(graphs4) or data(graphs5); see example. Currently, this package includes a list of graphs for three, four, or five lists.
#' @param logprior Log of the prior probability of each value in Nmissing. If left blank, this will default to the -log(Nmissing).
#' @param log.prior.model.weights Prior weights on the graphs. This should be a vector of the same length as `graphs`.
#' @param alpha Fractional posterior distribution hyperparameter. Defaults is 1 for regular posterior distirbutions.
#'
#' @return This function returns a matrix of weights, where rows correspond to models and columns correspond to values of Nmissing. Thus, the \code{ij}th entry of the matrix is the posterior probability of the \code{i}th model and the \code{j}th entry of Nmissing. Row sums return posterior probabilities by model.Column sums return posterior probabilities by value of  Nmissing.
#'
#' @details This is the main function in this package.  It performs capture-recapture (or multiple systems estimation) using Bayesian model averaging as outlined in Madigan and York (1997).
#'
#' @author Adapted by Olivier Binette \email{olivier.binette@gmail.com} from the dga::bma.cr function of James Johndrow \email{james.johndrow@gmail.com} and Kristian Lum \email{kl@hrdag.org}
#'
#' @examples
#' library(dga)
#' 
#' #### 5 list example from M & Y ##########
#' delta <- .5
#' Y <- c(0, 27, 37, 19, 4, 4, 1, 1, 97, 22, 37, 25, 2, 1, 3, 5, 83, 36, 34, 18, 3, 5, 0, 2, 30, 5, 23, 8, 0, 3, 0, 2)
#' Y <- array(Y, dim = c(2, 2, 2, 2, 2))
#' Nmissing <- 1:300
#' N <- Nmissing + sum(Y)
#' data(graphs5)
#' weights <- bma.cr(Y, Nmissing, delta, graphs5)
#'
#' ##### 3 list example from M & Y #######
#' Y <- c(0, 60, 49, 4, 247, 112, 142, 12)
#' Y <- array(Y, dim = c(2, 2, 2))
#'
#' delta <- 1
#' a <- 13.14
#' b <- 55.17
#'
#' Nmissing <- 1:300
#' N <- Nmissing + sum(Y)
#'
#' logprior <- N * log(b) - (N + a) * log(1 + b) + lgamma(N + a) - lgamma(N + 1) - lgamma(a)
#'
#' data(graphs3)
#' weights <- bma.cr(Y, Nmissing, delta, graphs3, logprior)
#'
#' @source `dga` package.
#' @references James Johndrow, Kristian Lum and Patrick Ball (2015). dga: Capture-Recapture Estimation using Bayesian Model Averaging. R package version 1.2. https://CRAN.R-project.org/package=dga
#' @export
bma.cr <- function(Y, Nmissing, delta, graphs,
                   logprior = NULL,
                   log.prior.model.weights = NULL,
                   alpha = 1) {
  UseMethod("bma.cr", Y)
}

#' @import assert
#' @export
bma.cr.MSEdata <- function(Y, Nmissing, delta, graphs,
                           logprior = NULL,
                           log.prior.model.weights = NULL,
                           alpha = 1) {
  assert(is.MSEdata(Y))

  Y <- MSEdata_to_array(Y)

  bma.cr.array(Y, Nmissing, delta, graphs,
    logprior = logprior,
    log.prior.model.weights = log.prior.model.weights,
    alpha = alpha
  )
}

#' @export
bma.cr.array <- function(Y, Nmissing, delta, graphs,
                         logprior = NULL,
                         log.prior.model.weights = NULL,
                         alpha = 1) {
  if (is.null(logprior)) {
    logprior <- -log(sum(Y) + Nmissing)
  }
  if (length(delta) == 1) {
    delta <- rep(delta, length(Y))
    delta <- array(delta, dim = dim(Y))
  }

  Y[1] <- 0
  p <- length(dim(Y))
  s_delta <- sum(delta)

  # Precomputations
  compMat <- MakeCompMatrix(p, delta, alpha * Y, alpha * Nmissing)
  D <- lgamma(s_delta) - lgamma(alpha * (Nmissing + sum(Y)) + s_delta)
  multinomialCoefficient <- alpha * (lgamma(Nmissing + sum(Y) + 1) - sum(lgamma(Y[-1] + 1)) - lgamma(Nmissing + 1))

  # Compute log posterior for all models
  weights <- computeLogPostProbs(compMat, graphs, D, p)
  rowAdd(weights, multinomialCoefficient)
  rowAdd(weights, logprior)
  if (!is.null(log.prior.model.weights)) colAdd(weights, log.prior.model.weights)

  # Normalization
  expNormalize(weights)

  return(weights)
}
