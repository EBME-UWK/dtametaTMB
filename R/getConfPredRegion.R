#' @keywords internal
#' @importFrom stats qlogis plogis qf
#' @noRd
getConfPredRegion <- function(muA,muB,
                              seA,seB,covAB, # confregion
                              varA,varB,sAB, # predregion
                              nstudy,
                              conflevel,
                              predlevel){  
  r       <- covAB / (seA*seB)
  sepredA <- sqrt(varA + seA**2)
  sepredB <- sqrt(varB + seB**2)
  rpredAB <- (sAB + covAB) / (sepredA*sepredB)
  f_conf  <- stats::qf(conflevel, df1 = 2, df2 = nstudy - 2)
  f_pred  <- stats::qf(predlevel, df1 = 2, df2 = nstudy - 2)
  croot_conf <- sqrt(2 * f_conf)
  croot_pred <- sqrt(2 * f_pred)
  
  conf_region <- c()
  pred_region <- c()
  # Confidence region
  for (i in seq(0, 2*pi, length.out=361)){
    confA    <- muA + (seA*croot_conf*cos(i))
    confB    <- muB + (seB*croot_conf*cos(i + acos(r)))
    confsens <- stats::plogis(confA)
    confspec <- stats::plogis(confB)
    conf_i   <- data.frame(X=1-confspec, Y=confsens)
    conf_region <- rbind(conf_region, conf_i)
  }
  for (i in seq(0, 2*pi, length.out=361)){
    predA    <- muA + (sepredA*croot_pred*cos(i))
    predB    <- muB + (sepredB*croot_pred*cos(i + acos(rpredAB)))
    predsens <- stats::plogis(predA)
    predspec <- stats::plogis(predB)
    pred_i   <- data.frame(X=1-predspec, Y=predsens)
    pred_region<-rbind(pred_region, pred_i)
  }
  region <- list(conf=conf_region,pred=pred_region)
  return(region)
}

