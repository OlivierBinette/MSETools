#' MSE data format
#'
#' The function \code{MSEdata()} transforms an existing dataframe to the "MSE" format,
#' ensuring it contains a "count" column and that the other columns refer to
#' inclusion (1) or exclusion (0) on a set of lists.
#'
#' Zero counts of unobserved capture patterns are added and duplicates capture patterns
#' are aggregated.
#'
#' @param data Original MSE dataframe. It should contain a column named representing
#'     the observed counts of capture patterns, as well as columns representing
#'     the different lists, as follows:
#' @param colname Name of the column representing the observed counts. Default is "count".
#' \preformatted{         c1    c2 count
#'          0     1     7
#'          1     0     3
#'          1     1     4}
#'
#' @import assert
#' @export
MSEdata <- function(data, colname = "count") {
  assert(inherits(data, "data.frame"))

  # Validate count column
  assert(colname %in% names(data),
         msg = paste0("A column named '", colname, "' should be specified.")
  )
  if (colname != "count") {
    colnames(data)[colnames(data) == colname] <- "count"
  }
  assert(is.numeric(data$count),
         all(data$count >= 0),
         all((data$count %% 1) == 0),
         msg = "Count column should only contain non-negative integers."
  )

  # Validate other columns
  listnames <- setdiff(names(data), "count")
  for (list in listnames) {
    assert(is.numeric(data[, list, drop = TRUE]))
    assert(all(data[, list, drop = TRUE] %in% c(0, 1)),
           msg = "List columns can only contain zeros and ones."
    )
  }

  data <- clean_MSE_data(data)

  attr(data, "class") <- c("MSEdata", attr(data, "class"))

  return(data)
}

#' Standardize MSE data format
#'
#' @param data MSE data frame.
#'
#' @importFrom dplyr left_join group_by_at vars count ungroup `%>%`
#' @importFrom purrr map_dfc
clean_MSE_data <- function(data) {
  nlists <- ncol(data) - 1
  data <- data %>%
    group_by_at(vars(-count)) %>%
    count(wt = count, name = "count") %>%
    ungroup()

  # Binary table with all combinations of zeros and ones
  X <- eval(parse(
    text =
      paste0("table(", paste0(rep("c(0,1)", nlists), collapse = ","), ")")
  )) %>%
    as.data.frame.table() %>%
    map_dfc(as.numeric) - 1

  # Removing the count for unobserved cases and removing superfluous column
  X <- X[2:nrow(X), 1:nlists]

  # Match column names of the data to those of the binary matrix
  listnames <- setdiff(names(data), "count")
  colnames(X) <- listnames

  # Join the binary table with the observed counts
  result <- left_join(X, data, by = listnames)

  # Reorder observations
  o1 <- order(rowApply(result, function(x) paste0(x, collapse = "")))
  result <- result[rev(o1), ]
  o2 <- order(rowSums(result[, listnames]))
  result <- result[o2, ]

  # Set NA counts to zero
  result[is.na(result[, "count"]), "count"] <- 0

  rownames(result) <- 1:nrow(result)

  return(result)
}

#' Check that data is of class `MSEdata`.
#'
#' @usage is.MSEdata(data)
#'
#' @param data Object.
#'
#' @export
is.MSEdata <- function(data) {
  inherits(data, "MSEdata")
}

#' Get list names
#'
#' Returns the names of the MSE lists.
#'
#' @usage listnames(data)
#'
#' @param data Data of type `MSEdata`.
#'
#' @export
listnames <- function(data) {
  assert(is.MSEdata(data))

  return(base::setdiff(names(data), "count"))
}

#' Set list names
#'
#' @usage `listnames<-`(data, value)
#'
#' @param data Data of type `MSEdata`.
#' @param value character vector of list names of length `nlists()`
#'
#' @export
`listnames<-` <- function(data, value) {
  assert(is.MSEdata(data))
  assert(length(value) == ncol(data) - 1)

  colnames(data)[colnames(data) != "count"] <- value
  return(data)
}

#' Number of lists
#'
#' This is ncol - 1.
#'
#' @usage nlists(data)
#'
#' @param data Data of type `MSEdata`.
#'
#' @export
nlists <- function(data) {
  assert(is.MSEdata(data))

  return(length(listnames(data)))
}


#' Omit lists from MSE data
#'
#' @usage omit(data, lists)
#'
#' @param data Data of type `MSEdata`.
#' @param lists Character vector of lists to be omitted.
#'
#' @export
omit <- function(data, lists) {
  assert(is.MSEdata(data))

  cols <- setdiff(names(data), lists)

  return(MSEdata(data[, cols]))
}

#' Merge lists
#'
#' @usage mergeLists(data, ...)
#'
#' @param data Data of type `MSEdata`.
#' @param ... Character vectors representing sets of lists to be merged.
#'
#' @export
mergeLists <- function(data, ...) {
  assert(is.MSEdata(data))

  args <- list(...)
  for (lists in args) {
    data[paste(lists, collapse = "-")] <- 1 * (rowSums(data[, lists]) > 0)
  }
  data <- data[, setdiff(names(x), unlist(args))]

  return(MSEdata(data))
}

MSEdata_to_array <- function(data) {
  Y <- data[, rev(which(colnames(data) != "count"))]
  o <- order(rowApply(Y, function(x) paste0(x, collapse = "")))
  Y <- array(c(0, data[o, "count"]), dim = rep(2, times = nlists(data)))
  return(Y)
}
