#' Relativize abundances by column
#'
#' This function takes a numeric matrix or data frame and divides each value
#' in a column by the column sum, so that each column sums to 1.
#' Useful for relativizing species abundances or similar data prior to performing
#' clustering or DISPROF analyses of species data.
#'
#' @param Y A numeric matrix or data frame.
#'
#' @return An object of the same class as `Y` with values scaled so that
#' each column sums to 1.
#'
#' @examples
#' Y <- data.frame(a = c(1, 2, 3), b = c(4, 5, 6))
#' relativize_cols(Y)
#'
#' @export
f_relativize_counts <- function(x) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("Input must be a matrix or data frame.")
  }
  if (!all(sapply(x, is.numeric))) {
    stop("All columns must be numeric.")
  }
  sweep(x, 2, colSums(x), "/")
}
