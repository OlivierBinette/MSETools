rowApply <- function(data, fun) {
  apply(data, 1, function(x) fun(as.vector(x)))
}
