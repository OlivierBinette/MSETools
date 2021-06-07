
#' @export
tex <- function(x, ...) {
  UseMethod("tex", x)
}

#' Export MSEdata object to Tex table
#'
#' @import assert dplyr
#' @export
tex.MSEdata <- function(x, file) {
  assert(is.MSEdata(x))

  x <- x %>% filter(count > 0)
  counts <- x$count
  x <- x %>%
    mutate_all(function(x) ifelse(x == 1, "$\\times$", "")) %>%
    t()
  colnames(x) <- counts

  table <- kableExtra::kbl(x, format = "latex", escape = FALSE, booktabs = TRUE)
  kableExtra::save_kable(table, file)
}
