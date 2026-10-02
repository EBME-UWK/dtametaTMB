#' @keywords internal
#' @noRd
getWEIGHTS <- function(xdata,size){
  if(size=="equal"){
    se <- rep(1,nrow(xdata))
    sp <- rep(1,nrow(xdata))
  }
  if(size=="sampsize"){
    se <- xdata$n1 / sum(xdata$n1)*100
    sp <- xdata$n0 / sum(xdata$n0)*100
  }
  if(size=="se"){
    sem1  <- xdata$sens*(1-xdata$sens)*xdata$n1 # inverse logit variance
    spm1  <- xdata$spec*(1-xdata$spec)*xdata$n0 # inverse logit variance
    sem1[is.nan(sem1)] <- 0
    spm1[is.nan(spm1)] <- 0
    se <- sqrt(sem1) / sum(sqrt(sem1))*100
    sp <- sqrt(spm1) / sum(sqrt(spm1))*100
  }
  if(size=="sampsize_revman"){
    se  <- (4 + 0.64*sqrt(xdata$n1))/540
    sp  <- (4 + 0.64*sqrt(xdata$n0))/540
  }
  if(size=="se_revman"){
    sem1  <- xdata$sens*(1-xdata$sens)*xdata$n1 # inverse logit variance
    spm1  <- xdata$spec*(1-xdata$spec)*xdata$n0 # inverse logit variance
    sem1[is.nan(sem1)] <- 0
    spm1[is.nan(spm1)] <- 0
    se  <- (4 + 1.5*sqrt(sem1))/540
    sp  <- (4 + 1.5*sqrt(spm1))/540
  }
  pct <- data.frame(sp=sp,
                    se=se)
  return(pct)
}
