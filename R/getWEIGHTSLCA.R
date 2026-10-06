#' Calculate graphical sizing measures for latent-class SROC plots
#'
#' Internal helper calculating horizontal specificity and vertical
#' sensitivity dimension measures for study symbols in latent-class
#' SROC plots.
#'
#' The empirical Bayes methods use inverse conditional posterior standard
#' errors of the study-specific logit-sensitivity and logit-specificity
#' estimates. Non-finite or negative inverse posterior variances are assigned
#' zero for graphical sizing.
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
#' For \code{size = "eb_revman"}, the raw measures are inverse conditional
#' posterior standard errors of the empirical Bayes logit-specificity and
#' logit-sensitivity estimates. This is an analogous graphical extension
#' for latent-class models and is not an original RevMan latent-class
#' weighting method.
#'
#' For \code{size = "sampsize_revman"}, both graphical dimensions are based
#' on the square root of the total study sample size. Consequently, this
#' method produces symbols with equal horizontal and vertical dimensions.
#' A total sample size of zero is replaced by one for graphical sizing.
#'
#' Unlike the cell-count-based RevMan calculations used for conventional
#' DTA models, these latent-class methods do not apply zero-cell continuity
#' corrections. All calculations and adjustments affect graphical dimensions
#' only and do not affect model fitting.
#'
#' @param xdata A data frame containing total study sample sizes and
#'   conditional posterior variances of the empirical Bayes
#'   logit-sensitivity and logit-specificity estimates.
#' @param size A character string identifying the symbol-sizing method.
#'
#' @return A data frame with columns \code{sp} and \code{se}, containing
#'   horizontal specificity and vertical sensitivity graphical dimensions,
#'   respectively.
#'
#' @keywords internal
#' @noRd
getWEIGHTSLCA <- function(xdata,size){
  if(size=="eb"){
    invvar_sens <- 1/xdata$lsens_eb_var # inverse logit eb variance
    invvar_spec <- 1/xdata$lspec_eb_var # inverse logit eb variance
    invvar_sens[!is.finite(invvar_sens) | invvar_sens < 0] <- 0
    invvar_spec[!is.finite(invvar_spec) | invvar_spec < 0] <- 0
    invse_sens <- sqrt(invvar_sens)
    invse_spec <- sqrt(invvar_spec)
    if (sum(invse_sens) > 0) { se <- invse_sens / sum(invse_sens) * 100 
     } else { se <- rep(1, nrow(xdata))}
    if (sum(invse_spec) > 0) { sp <- invse_spec / sum(invse_spec) * 100
     } else { sp <- rep(1, nrow(xdata))}
  }
  if(size=="equal"){
    se <- rep(1,nrow(xdata))
    sp <- rep(1,nrow(xdata))
  }
  if(size=="sampsize"){
    se <- xdata$n / sum(xdata$n)*100
    sp <- xdata$n / sum(xdata$n)*100
  }
  if(size=="eb_revman") {
    invvar_sens <- 1/xdata$lsens_eb_var # inverse logit eb variance
    invvar_spec <- 1/xdata$lspec_eb_var # inverse logit eb variance
    invvar_sens[!is.finite(invvar_sens) | invvar_sens < 0] <- 0
    invvar_spec[!is.finite(invvar_spec) | invvar_spec < 0] <- 0
    invse_sens <- sqrt(invvar_sens)
    invse_spec <- sqrt(invvar_spec)
    c_eb   <- 8/mean(pmax(invse_sens,invse_spec))
    if(!is.finite(c_eb)){c_eb <- 0}
    se  <- (4 + c_eb*invse_sens)/10 # revman_rescale: 10
    sp  <- (4 + c_eb*invse_spec)/10 # revman_rescale: 10
  }
  if(size=="sampsize_revman"){
    n       <- xdata$n
    n[n==0] <- 1
    sn      <- sqrt(n)
    c_sn    <- 8/mean(sn)
    se   <- (4 + c_sn*sn)/10 # revman_rescale: 10
    sp   <- (4 + c_sn*sn)/10 # revman_rescale: 10
  }
  pct <- data.frame(sp=sp,
                    se=se)
  return(pct)
}