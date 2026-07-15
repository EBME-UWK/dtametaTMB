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
  pct <- data.frame(sp=pctsp,
                    se=pctse)
  return(pct)
}