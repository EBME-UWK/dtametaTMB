#' @keywords internal
#' @noRd
getXP <- function(X){
  X$n1    <- X$TP+X$FN
  X$n0    <- X$FP+X$TN
  X$sens  <- X$TP / X$n1
  X$spec  <- X$TN / X$n0
  X$recordid <- seq_len(nrow(X))
  return(X)
}
