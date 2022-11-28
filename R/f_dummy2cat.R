#' Convert dummy codes to a single categorical variable \cr
#' \cr
#' The purpose of this code is to essentially undo the effects of f_dummy
#'
#' @param X binary dummy codes created by f_dummy
#' @return column vector of integers specifying factor levels or group
#'   membership

f_dummy2cat <- function(X){
  X_nR <- dim(X)[1]
  X_nC <- dim(X)[2]
  temp <- matrix(0, X_nR, 1)
  for (i in 1:X_nC) {
    temp[X[,i]==TRUE] <- i
  }
  return(temp)
}
