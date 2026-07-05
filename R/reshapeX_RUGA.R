#' @keywords internal
#' @importFrom stats reshape
#' @noRd
reshapeX_RUGA <- function(X){
  X$y1 <- X$TP
  X$y0 <- X$FP  
  X$n1 <- X$TP+X$FN
  X$n0 <- X$FP+X$TN
  X$recordid <- seq_len(nrow(X)) 
  ### Reshape the data from wide to long format. ###
  Y = stats::reshape(X, direction="long", varying=list(c("n1", "n0"), c("y1","y0")),
                     timevar="x", times=c(0.5,-0.5), v.names=c("n","y")) 
  Y = Y[order(Y$recordid),] 
  return(Y)
}
