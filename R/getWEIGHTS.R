#' Calculate graphical sizing measures for SROC study symbols
#'
#' Internal helper calculating horizontal specificity and vertical
#' sensitivity dimension measures for the supported symbol-sizing methods.
#'
#' The RevMan-style methods use full dimensions of the form
#'
#' d_i = 4 + c * q_i,
#'
#' where c is chosen so that the mean study-level maximum of c * q is 8.
#' The resulting dimensions are divided by 10 to place them on the package's
#' graphical scale. No continuity correction is applied to the inverse-SE
#' measures.
#'
#' @param xdata Data frame containing study-level sensitivity, specificity,
#' and diseased and non-diseased sample sizes.
#' @param size Character string identifying the symbol-sizing method.
#'
#' @return A data frame with \code{sp} and \code{se}, containing horizontal
#' and vertical graphical dimension measures, respectively.
#'
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
    q_se <- sqrt(xdata$n1)
    q_sp <- sqrt(xdata$n0)
    c_ss <- 8 / mean(pmax(q_se,q_sp))
    se  <- (4 + c_ss*q_se)/10 # revman_scale: 10
    sp  <- (4 + c_ss*q_sp)/10 # revman_scale: 10
  }
  if(size=="se_revman"){
    sem1  <- xdata$sens*(1-xdata$sens)*xdata$n1 # inverse logit variance
    spm1  <- xdata$spec*(1-xdata$spec)*xdata$n0 # inverse logit variance
    sem1[is.nan(sem1)] <- 0
    spm1[is.nan(spm1)] <- 0
    q_se  <- sqrt(sem1)
    q_sp  <- sqrt(spm1)
    c_se  <- 8/mean(pmax(q_se,q_sp))
    se  <- (4 + c_se*q_se)/10 # revman_scale: 10
    sp  <- (4 + c_se*q_sp)/10 # revman_scale: 10
  }
  pct <- data.frame(sp=sp,
                    se=se)
  return(pct)
}