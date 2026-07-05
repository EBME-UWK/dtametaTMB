#' @keywords internal
#' @noRd
reshapeX_REIT <- function(X){
  X$true1 <- X$TP
  X$true0 <- X$TN 
  X$n1    <- X$TP+X$FN
  X$n0    <- X$FP+X$TN
  X$recordid <- seq_len(nrow(X))
  Y <- reshape(X, direction="long", varying=list(c("n1", "n0"), c("true1", "true0")), 
               timevar="sens", times=c(1,0), v.names=c("n","true")) 
  Y <- Y[order(Y$recordid),]  
  Y$spec <- 1-Y$sens
  return(Y)
}
