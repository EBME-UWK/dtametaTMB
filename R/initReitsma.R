#' @keywords internal
#' @noRd

initReitsma <- function(X) {
  has_sens_info   <- with(X,TP + FN > 0)
  has_spec_info   <- with(X,TN + FP > 0)
  has_paired_info <- has_sens_info & has_spec_info
  logit_sens   <- with(X[has_sens_info,,drop=FALSE],log((TP+0.5)/(FN+0.5)))
  logit_spec   <- with(X[has_spec_info,,drop=FALSE],log((TN+0.5)/(FP+0.5)))
  ##
  muA_init     <- mean(logit_sens,na.rm=TRUE)
  if(!is.finite(muA_init)) { muA_init <- 0 }
  muB_init     <- mean(logit_spec,na.rm=TRUE)
  if(!is.finite(muB_init)) { muB_init <- 0 }
  ##
  sA_init      <- stats::sd(logit_sens,na.rm=TRUE)
  if (!is.finite(sA_init)) {sA_init <- 0.5}
  sA_init      <- max(sA_init,1e-05)
  ##
  sB_init      <- stats::sd(logit_spec,na.rm=TRUE)
  if (!is.finite(sB_init)) {sB_init <- 0.5}
  sB_init      <- max(sB_init,1e-05)
  ##
  logit_sens_cor <- with(X[has_paired_info,,drop=FALSE],log((TP+0.5)/(FN+0.5)))
  logit_spec_cor <- with(X[has_paired_info,,drop=FALSE],log((TN+0.5)/(FP+0.5)))
  rAB_init <- 0
  if (length(logit_sens_cor) >= 3L) {
    sd_sens_cor <- stats::sd(logit_sens_cor)
    sd_spec_cor <- stats::sd(logit_spec_cor)
    if (is.finite(sd_sens_cor) && is.finite(sd_spec_cor) && sd_sens_cor > sqrt(.Machine$double.eps) && sd_spec_cor > sqrt(.Machine$double.eps)) {
      empirical_cor <- suppressWarnings(stats::cor(logit_sens_cor,logit_spec_cor))
      if (is.finite(empirical_cor)) {rAB_init <- min(max(empirical_cor, -0.99),0.99)}
    }
  }
  return(list(muA_init=muA_init,
              muB_init=muB_init,
              sA_init=sA_init,
              sB_init=sB_init,
              rAB_init=rAB_init))
}

