#' Recode elements of vector as consecutive integers \cr
#' \cr
#' This function is used to convert treatment levels (e.g., "[2 2 2 5 5 5]") to
#' "[1 1 1 2 2 2]" so proper ANOVA design matrices can be constructed.
#'
#' @param Y vector or matrix of input variables (column-wise)
#' @return vector or matrix recoded as consecutive integers

f_recode <- function(Y){
  Y_nR <- length(Y)
  X    <- rep(NaN, Y_nR)   # Preallocate matrix
  G    <- sort(unique(Y))  # Get sorted list of unique values
  G_nR <- length(G)        # Get no. of unique values
  temp <- rep(NaN, G_nR)   # Preallocate vector
  for (j in 1:G_nR) {
    X[Y==G[j]] <- j
  }
  return(X)
}
