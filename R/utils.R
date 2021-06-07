
#' Compute adjacency matrix of a decomposable graph
#'
#' @param graph Decomposable graph as in `data(graphs5)`.
#' @param p of vertices of the graph.
#'
#' @return Adjacency matrix.
#'
adjMat <- function(graph, p) {
  mat <- matrix(rep(0, p^2), nrow = p, ncol = p)
  for (clique in graph$C) {
    clique <- c(clique)
    mat[clique, clique] <- 1
  }
  return(mat - diag(1, p))
}

CompLogML <- function(D, Nmissing, delta) {
  Nmissing <- Nmissing + D[1] + delta[1]
  lgamma(Nmissing) + sum(lgamma(D[2:length(D)] + delta[2:length(D)])) - sum(lgamma(delta))
}

integer.base.b <- function(x, b = 2) {
  xi <- as.integer(x)
  if (any(is.na(xi) | ((x - xi) != 0))) {
    print(list(ERROR = "x not integer", x = x))
  }
  N <- length(x)
  xMax <- max(x)
  ndigits <- (floor(logb(xMax, base = 2)) + 1)
  Base.b <- array(NA, dim = c(N, ndigits))
  for (i in 1:ndigits) {
    Base.b[, ndigits - i + 1] <- (x %% b)
    x <- (x %/% b)
  }
  if (N == 1) Base.b[1, ] else Base.b
}

MakeCompMatrix <- function(p, delta, Y, Nmissing) {
  compLMLs <- matrix(0, nrow = 2^p - 1, ncol = length(Nmissing))
  bins <- integer.base.b(1:(2^p - 1), 2)
  for (i in 1:(2^p - 1)) {
    inds <- which(bins[i, ] == 1)
    D <- c(apply(Y, inds, sum))
    alpha <- c(apply(delta, inds, sum))
    compLMLs[i, ] <- CompLogML(D, Nmissing, alpha)
  }
  return(compLMLs)
}

rowApply <- function(data, fun) {
  apply(data, 1, function(x) fun(as.vector(x)))
}
