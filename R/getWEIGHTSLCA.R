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
    pctse  <- (4 + 1.5*invse_sens)/5 # revman_rescale: 5
    pctsp  <- (4 + 1.5*invse_spec)/5 # revman_rescale: 5
  }
  if(size=="sampsize_revman"){
    pctse  <- (4 + 0.64*sqrt(xdata$n))/5 # revman_rescale: 5
    pctsp  <- (4 + 0.64*sqrt(xdata$n))/5 # revman_rescale: 5
  }
  pct <- data.frame(sp=pctsp,
                    se=pctse)
  return(pct)
}