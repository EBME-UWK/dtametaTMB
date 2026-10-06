#' Calculate graphical sizing measures for SROC study symbols
#'
#' Internal helper calculating horizontal specificity and vertical
#' sensitivity dimension measures for the supported symbol-sizing methods.
#'
#' For the RevMan-style methods, the raw horizontal and vertical measures
#' are converted to full graphical dimensions according to
#'
#' \deqn{
#' d_{ij} = \frac{4 + c q_{ij}}{10},
#' }
#'
#' where a common multiplier \eqn{c} is chosen such that the mean
#' study-level maximum of \eqn{c q_{ij}} is 8:
#'
#' \deqn{
#' c =
#' \frac{8}{
#' \operatorname{mean}_i
#' \left\{
#' \max(q_{i,\mathrm{sp}}, q_{i,\mathrm{se}})
#' \right\}
#' }.
#' }
#'
#' For \code{size = "sampsize_revman"}, the raw dimensions are the square
#' roots of the non-diseased and diseased sample sizes. If either sample
#' size is zero, 1 is added to both sample sizes before calculating the raw
#' dimensions, following the RevMan calculation.
#'
#' For \code{size = "se_revman"}, the raw dimensions are approximate inverse standard
#' errors of the empirical logit specificity and sensitivity. If any of the
#' four cell counts is zero, 0.5 is added to all four cells before calculating
#' the raw dimensions, following the RevMan calculation.
#'
#' No continuity correction is applied for \code{size = "se"}. All
#' adjustments described above affect graphical dimensions only and do not
#' alter the cell counts used for model fitting.
#'
#' @param xdata A data frame containing study-level cell counts,
#' sensitivities, specificities, and diseased and non-diseased sample
#' sizes.
#' @param size A character string identifying the symbol-sizing method.
#'
#' @return A data frame with columns \code{sp} and \code{se}, containing
#' horizontal specificity and vertical sensitivity graphical dimensions,
#' respectively.
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
    sem1  <- xdata$sens*(1-xdata$sens)*xdata$n1 # inverse variance of empirical logit sensitivity
    spm1  <- xdata$spec*(1-xdata$spec)*xdata$n0 # inverse variance of empirical logit specificity
    sem1[!is.finite(sem1)] <- 0
    spm1[!is.finite(spm1)] <- 0
    sqrtsem1 <- sqrt(sem1)
    sqrtspm1 <- sqrt(spm1)
    if (sum(sqrtsem1) > 0) { se <- sqrtsem1 / sum(sqrtsem1) * 100 
     } else { se <- rep(1, nrow(xdata))}
    if (sum(sqrtspm1) > 0) { sp <- sqrtspm1 / sum(sqrtspm1) * 100
     } else { sp <- rep(1, nrow(xdata))}
  }
  if(size=="sampsize_revman"){
    n1 <- xdata$n1
    n0 <- xdata$n0
    zero_dimension <- ( n1 == 0 | n0 == 0 )
    n1[zero_dimension] <- n1[zero_dimension] + 1
    n0[zero_dimension] <- n0[zero_dimension] + 1
    q_se <- sqrt(n1)
    q_sp <- sqrt(n0)
    c_ss <- 8 / mean(pmax(q_se, q_sp))
    se <- ( 4 + c_ss * q_se ) / 10 # revman_scale: 10
    sp <- ( 4 + c_ss * q_sp ) / 10 # revman_scale: 10
  }
  if(size=="se_revman") {
    TP <- xdata$TP
    FP <- xdata$FP
    FN <- xdata$FN
    TN <- xdata$TN
    zero_cell <- (TP == 0 | FP == 0 | FN == 0 | TN == 0)
    TP[zero_cell] <- TP[zero_cell] + 0.5
    FP[zero_cell] <- FP[zero_cell] + 0.5
    FN[zero_cell] <- FN[zero_cell] + 0.5
    TN[zero_cell] <- TN[zero_cell] + 0.5
    q_se <- 1 / sqrt( 1 / TP + 1 / FN)
    q_sp <- 1 / sqrt( 1 / TN + 1 / FP)
    c_se <- 8 /  mean(pmax(q_se, q_sp))
    se <- ( 4 + c_se * q_se ) / 10 # revman_scale: 10
    sp <- ( 4 + c_se * q_sp ) / 10 # revman_scale: 10
  }
  pct <- data.frame(sp=sp,
                    se=se)
  return(pct)
}