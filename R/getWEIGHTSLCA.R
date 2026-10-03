#' Calculate graphical sizing measures for latent-class SROC plots
#'
#' Internal helper calculating horizontal specificity and vertical
#' sensitivity dimension measures for study symbols in latent-class
#' SROC plots.
#'
#' The empirical Bayes methods use inverse conditional posterior standard
#' errors of the study-specific logit-sensitivity and logit-specificity
#' estimates. The RevMan-style methods use full dimensions of the form
#'
#' d_i = 4 + c * q_i,
#'
#' where c is chosen so that the mean study-level maximum of c * q is 8.
#' The resulting dimensions are divided by 10 to place them on the
#' package's graphical scale.
#'
#' The empirical Bayes RevMan-style method is an analogous graphical
#' extension for latent-class models and is not an original RevMan
#' latent-class weighting method.
#'
#' @param xdata Data frame containing total study sample size and conditional
#' posterior variances of the empirical Bayes logit-sensitivity and
#' logit-specificity estimates.
#' @param size Character string identifying the symbol-sizing method.
#'
#' @return A data frame with \code{sp} and \code{se}, containing horizontal
#' and vertical graphical dimension measures, respectively.
#'
#' @keywords internal
#' @noRd
getWEIGHTSLCA <- function(xdata,size){
  if(size=="eb"){
    se_vu <- 1/xdata$lsens_eb_var # inverse logit eb variance
    sp_vu <- 1/xdata$lspec_eb_var # inverse logit eb variance
    pctse <- sqrt(se_vu) / sum(sqrt(se_vu))*100
    pctsp <- sqrt(sp_vu) / sum(sqrt(sp_vu))*100
  }
  if(size=="equal"){
    pctse <- rep(1,nrow(xdata))
    pctsp <- rep(1,nrow(xdata))
  }
  if(size=="sampsize"){
    pctse <- xdata$n / sum(xdata$n)*100
    pctsp <- xdata$n / sum(xdata$n)*100
  }
  if(size=="eb_revman") {
    invse_sens <- 1 / sqrt(xdata$lsens_eb_var)
    invse_spec <- 1 / sqrt(xdata$lspec_eb_var)
    c_eb   <- 8/mean(pmax(invse_sens,invse_spec))
    pctse  <- (4 + c_eb*invse_sens)/10 # revman_rescale: 10
    pctsp  <- (4 + c_eb*invse_spec)/10 # revman_rescale: 10
  }
  if(size=="sampsize_revman"){
    sn     <- sqrt(xdata$n)
    c_sn   <- 8/mean(sn)
    pctse  <- (4 + c_sn*sn)/10 # revman_rescale: 10
    pctsp  <- (4 + c_sn*sn)/10 # revman_rescale: 10
  }
  pct <- data.frame(sp=pctsp,
                    se=pctse)
  return(pct)
}