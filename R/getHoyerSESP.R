#' @keywords internal
#' @noRd
getHoyerSESP <- function(rep2, threshold, testdir, conflevel) {
  qq <- qnorm(1-(1-conflevel)/2)
  rls1 <- which(rownames(rep2)=="logitSurv1")
  rls0 <- which(rownames(rep2)=="logitSurv0")
  
  if(testdir=="greater"){
    sens <- data.frame(threshold=threshold,
                       logitSurv1=rep2[rls1,"Estimate"],
                       Std_Error=rep2[rls1,"Std. Error"],
                       CI_Lower=NA,
                       CI_Upper=NA)
    sens$CI_Lower   <- with(sens,logitSurv1-qq*Std_Error)
    sens$CI_Upper   <- with(sens,logitSurv1+qq*Std_Error)
    sens$Sens       <- with(sens,plogis(logitSurv1))
    sens$SensCI_Lower <- with(sens,plogis(CI_Lower))
    sens$SensCI_Upper <- with(sens,plogis(CI_Upper))
    
    spec <- data.frame(threshold=threshold,
                       logitSurv0=rep2[rls0,"Estimate"],
                       Std_Error=rep2[rls0,"Std. Error"],
                       CI_Lower=NA,
                       CI_Upper=NA)
    spec$CI_Lower   <- with(spec,logitSurv0-qq*Std_Error)
    spec$CI_Upper   <- with(spec,logitSurv0+qq*Std_Error)
    spec$Spec       <- with(spec,1-plogis(logitSurv0))
    spec$SpecCI_Lower <- with(spec,1-plogis(CI_Upper))
    spec$SpecCI_Upper <- with(spec,1-plogis(CI_Lower))
    
  }
  ########
  if(testdir=="less"){
    sens <- data.frame(threshold=threshold,
                       logitSurv1=rep2[rls1,"Estimate"],
                       Std_Error=rep2[rls1,"Std. Error"],
                       CI_Lower=NA,
                       CI_Upper=NA)
    sens$CI_Lower   <- with(sens,logitSurv1-qq*Std_Error)
    sens$CI_Upper   <- with(sens,logitSurv1+qq*Std_Error)
    sens$Sens       <- with(sens,1-plogis(logitSurv1))
    sens$SensCI_Lower <- with(sens,1-plogis(CI_Upper))
    sens$SensCI_Upper <- with(sens,1-plogis(CI_Lower))
    
    spec <- data.frame(threshold=threshold,
                       logitSurv0=rep2[rls0,"Estimate"],
                       Std_Error=rep2[rls0,"Std. Error"],
                       CI_Lower=NA,
                       CI_Upper=NA)
    spec$CI_Lower   <- with(spec,logitSurv0-qq*Std_Error)
    spec$CI_Upper   <- with(spec,logitSurv0+qq*Std_Error)
    spec$Spec       <- with(spec,plogis(logitSurv0))
    spec$SpecCI_Lower <- with(spec,plogis(CI_Lower))
    spec$SpecCI_Upper <- with(spec,plogis(CI_Upper))
    
  }
  
  sesp <- data.frame(threshold=sens$threshold,
                     conflevel=conflevel,
                     Sens=sens$Sens,
                     SensCI_Lower=sens$SensCI_Lower,
                     SensCI_Upper=sens$SensCI_Upper,
                     Spec=spec$Spec,
                     SpecCI_Lower=spec$SpecCI_Lower,
                     SpecCI_Upper=spec$SpecCI_Upper)
  return(sesp)
}