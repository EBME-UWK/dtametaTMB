#' @keywords internal
#' @importFrom stats qlogis plogis
#' @noRd
getROCpoints <- function(Lambda,beta,specrange){
  roc_points <- c()
  fromL  <- specrange[1]
  toU    <- specrange[2]
  for (i in seq(from=fromL, to=toU, by=0.005)){
    Sp_i  <- i
    Fpr_i <- 1-Sp_i
    LSp_i <- stats::qlogis(Sp_i)
    LSe_i <- Lambda*exp(-beta/2) - exp(-beta)*LSp_i
    Se_i  <- stats::plogis(LSe_i)
    roc_i <- data.frame(FPR=Fpr_i, Sen=Se_i)
    roc_points<-rbind(roc_points,roc_i)
  }
  return(roc_points)
}
