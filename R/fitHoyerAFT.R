#' Fit Threshold-Based bivariate time-to-event DTA Models
#'
#' Fits a bivariate accelerated failure time (AFT) model for diagnostic
#' test accuracy (DTA) data using interval-censored likelihoods as
#' described in Hoyer et al. (2018). The model is estimated using
#' Template Model Builder (TMB) with study-specific random effects.
#' 
#' The model reports logit-survival quantities, which correspond to sensitivity 
#' and false positive rate for \code{testdirection = "greater"}, and to false 
#' negative rate and specificity for \code{testdirection = "less"}.
#'
#' @param data A list as produced by \code{\link{restructure}}, containing:
#'   \describe{
#'     \item{restructured}{Interval-formatted data}
#'     \item{original}{Processed original data with derived quantities}
#'   }
#'
#' @param init A data frame of initial parameter values as produced by
#'   \code{\link{getInitParms}}.
#'
#' @param conflevel Confidence level for confidence intervals for sensitivities
#'   and specificities at the chosen thresholds. Defaults to \code{0.95}.
#'
#' @param threshold Optional numeric value or vector specifying the
#'   threshold(s) at which sensitivity and specificity should be evaluated.
#'   If \code{NA} (default), the median threshold from the original data
#'   is used.
#'
#' @return An object of class \code{"HoyerAFT"} containing:
#' \describe{
#'   \item{data}{Processed original data}
#'   \item{restructured}{Interval-formatted data}
#'   \item{fit}{Optimization output from \code{nlminb}}
#'   \item{sdreport}{TMB report object from \code{sdreport}}
#'   \item{sdreport2}{Summary of reported parameters}
#'   \item{distcode}{Distribution code used in the model}
#'   \item{sensspec}{Sensitivities and specificities at the chosen thresholds}
#' }
#'
#' @examples
#' \dontrun{
#' data("diabetes")
#' res <- restructure(
#'   data = diabetes,
#'   TP = TP,
#'   FP = FP,
#'   FN = FN,
#'   TN = TN,
#'   threshold = threshold,
#'   study = study,
#'   smallest = 2,
#'   largest = 10
#' )
#'
#' init <- getInitParms(res$restructured)
#'
#' fit <- fitHoyerAFT(res, init)
#'
#' summary(fit)
#' }
#'
#' @details
#' Model fitting is performed using \code{nlminb}, and uncertainty
#' estimates are obtained via \code{TMB::sdreport}. Parameter
#' transformations include log-scale parameters and Fisher's Z-transform
#' for correlations.
#' 
#' @references
#' Hoyer, A., Hirt, S., Kuss, O. (2018).
#' Meta-analysis of full ROC curves using bivariate time-to-event models for interval-censored data.
#' \emph{Research Synthesis Methods}, 9(1), 62-72.
#' \doi{10.1002/jrsm.1273}
#'
#' @importFrom TMB MakeADFun sdreport
#' @importFrom stats nlminb median qnorm
#'
#' @note Requires a compiled TMB model named \code{"Hoyer"}.
#' @export
fitHoyerAFT <- function(data, init, conflevel=0.95, threshold = NA) {

  # Extract components
  datao   <- data$original
  testdir <- unique(datao$testdirection)
  if (length(testdir) != 1) stop("testdirection must be unique")
  datar   <- data$restructured

  # Validate init
  required <- c("beta0_init", "lambda0_init",
                "beta1_init", "lambda1_init",
                "su0_init", "su1_init",
                "coru0u1_init", "distcode")

  if (!is.data.frame(init) || any(!required %in% names(init))) {
    stop("'init' must be a valid output from getInitParms().")
  }

  if (!all(is.na(threshold))) {

    if (!is.numeric(threshold)) {
      stop("'threshold' must be numeric.")
    }

    if (any(!is.finite(threshold))) {
      stop("'threshold' must contain only finite values.")
    }

    if (any(threshold <= 0)) {
      stop("'threshold' values must be positive.")
    }
  }

  # Prepare TMB data
  dat2 <- list(
    lowerB    = datar$lowerB,
    upperB    = datar$upperB,
    events0   = datar$events0,
    events1   = datar$events1,
    ctype     = as.integer(datar$ctype),
    threshold = ifelse(is.na(threshold),
                       stats::median(datao$threshold),
                       threshold),
    study     = as.integer(factor(datar$study)) - 1,
    nstudy    = length(unique(datar$study)),
    dist      = init$distcode
  )

  # Parameters
  parameters <- list(
    beta0       = init$beta0_init,
    log_lambda0 = log(init$lambda0_init),
    beta1       = init$beta1_init,
    log_lambda1 = log(init$lambda1_init),
    log_su0     = log(init$su0_init),
    log_su1     = log(init$su1_init),
    rho_trans   = atanh(init$coru0u1_init),
    u0 = rep(0, dat2$nstudy),
    u1 = rep(0, dat2$nstudy)
  )
  
  # TMB objective
  dat2$model <- "Hoyer"
  obj <- TMB::MakeADFun(data=dat2,
                        parameters,
                        random = c("u0", "u1"),
                        DLL = "dtametaTMB_TMBExports")

  # Optimization
  fit <- stats::nlminb(
    obj$par,
    obj$fn,
    obj$gr
  )

  # Reports
  rep  <- TMB::sdreport(obj)
  rep2 <- summary(rep, select = "report")
  
  # Get sensitivities and specificities at the thresholds
  qq <- qnorm(1-(1-conflevel)/2)
  rls1 <- which(rownames(rep2)=="logitSurv1")
  rls0 <- which(rownames(rep2)=="logitSurv0")
  
  if(testdir=="greater"){
    rep6 <- data.frame(threshold=dat2$threshold,
                       logitSurv1=rep2[rls1,"Estimate"],
                       Std_Error=rep2[rls1,"Std. Error"],
                       CI_Lower=NA,
                       CI_Upper=NA)
    rep6$CI_Lower   <- with(rep6,logitSurv1-qq*Std_Error)
    rep6$CI_Upper   <- with(rep6,logitSurv1+qq*Std_Error)
    rep6$Sens       <- with(rep6,plogis(logitSurv1))
    rep6$SensCI_Lower <- with(rep6,plogis(CI_Lower))
    rep6$SensCI_Upper <- with(rep6,plogis(CI_Upper))
  
    rep7 <- data.frame(threshold=dat2$threshold,
                       logitSurv0=rep2[rls0,"Estimate"],
                       Std_Error=rep2[rls0,"Std. Error"],
                       CI_Lower=NA,
                       CI_Upper=NA)
    rep7$CI_Lower   <- with(rep7,logitSurv0-qq*Std_Error)
    rep7$CI_Upper   <- with(rep7,logitSurv0+qq*Std_Error)
    rep7$Spec       <- with(rep7,1-plogis(logitSurv0))
    rep7$SpecCI_Lower <- with(rep7,1-plogis(CI_Upper))
    rep7$SpecCI_Upper <- with(rep7,1-plogis(CI_Lower))
   
    rep8 <- data.frame(threshold=rep6$threshold,
                       conflevel=conflevel,
                       Sens=rep6$Sens,
                       SensCI_Lower=rep6$SensCI_Lower,
                       SensCI_Upper=rep6$SensCI_Upper,
                       Spec=rep7$Spec,
                       SpecCI_Lower=rep7$SpecCI_Lower,
                       SpecCI_Upper=rep7$SpecCI_Upper)
  }
  ########
  if(testdir=="less"){
    rep6 <- data.frame(threshold=dat2$threshold,
                       logitSurv1=rep2[rls1,"Estimate"],
                       Std_Error=rep2[rls1,"Std. Error"],
                       CI_Lower=NA,
                       CI_Upper=NA)
    rep6$CI_Lower   <- with(rep6,logitSurv1-qq*Std_Error)
    rep6$CI_Upper   <- with(rep6,logitSurv1+qq*Std_Error)
    rep6$Sens       <- with(rep6,1-plogis(logitSurv1))
    rep6$SensCI_Lower <- with(rep6,1-plogis(CI_Upper))
    rep6$SensCI_Upper <- with(rep6,1-plogis(CI_Lower))
    
    rep7 <- data.frame(threshold=dat2$threshold,
                       logitSurv0=rep2[rls0,"Estimate"],
                       Std_Error=rep2[rls0,"Std. Error"],
                       CI_Lower=NA,
                       CI_Upper=NA)
    rep7$CI_Lower   <- with(rep7,logitSurv0-qq*Std_Error)
    rep7$CI_Upper   <- with(rep7,logitSurv0+qq*Std_Error)
    rep7$Spec       <- with(rep7,plogis(logitSurv0))
    rep7$SpecCI_Lower <- with(rep7,plogis(CI_Lower))
    rep7$SpecCI_Upper <- with(rep7,plogis(CI_Upper))
    
    rep8 <- data.frame(threshold=rep6$threshold,
                       conflevel=conflevel,
                       Sens=rep6$Sens,
                       SensCI_Lower=rep6$SensCI_Lower,
                       SensCI_Upper=rep6$SensCI_Upper,
                       Spec=rep7$Spec,
                       SpecCI_Lower=rep7$SpecCI_Lower,
                       SpecCI_Upper=rep7$SpecCI_Upper)
  }
  # Result object
  res <- list(
    data         = datao,
    restructured = datar,
    fit          = fit,
    sdreport     = rep,
    sdreport2    = rep2,
    distcode     = dat2$dist,
    sensspec     = rep8
  )

  # Assign class
  class(res) <- c("HoyerAFT","Cochrane")

  return(res)
}
