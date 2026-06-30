#' Fit Threshold-Based Bivariate Time-to-Event Model (Hoyer AFT)
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
#' @param data A list as produced by \code{\link{restructure_data}}, containing:
#'   \describe{
#'     \item{restructured}{Interval-formatted data}
#'     \item{original}{Processed original data with derived quantities}
#'   }
#'
#' @param init A data frame of initial parameter values as produced by
#'   \code{\link{initHoyerAFT}}.
#'
#' @param conflevel Confidence level for confidence intervals for sensitivities
#'   and specificities at the chosen thresholds. Defaults to \code{0.95}.
#'
#' @param threshold Optional numeric value or vector specifying the 
#'   prediction grid threshold(s) at which sensitivity and specificity 
#'   should be evaluated. If \code{NA} (default), the median threshold from 
#'   the original data is used.
#'   
#' @param verbose Whether TMB optimization output should be printed (default: FALSE).
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
#' data("diabetes")
#' res <- restructure_data(
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
#' init <- initHoyerAFT(res$restructured)
#' fit <- fitHoyerAFT(res, init)
#' summary(fit)
#' 
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
fitHoyerAFT <- function(data, init, conflevel=0.95, threshold = NA, verbose=FALSE) {

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
    stop("'init' must be a valid output from initHoyerAFT().")
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
  
  dat2$model <- "Hoyer"

  # TMB objective
  obj <- TMB::MakeADFun(data=dat2,
                        parameters,
                        random = c("u0", "u1"),
                        silent = !verbose,
                        DLL = "dtametaTMB_TMBExports")
  # Optimization
  fit <- stats::nlminb(obj$par,
                       obj$fn,
                       obj$gr)
  
  # Convergence warning.
  if (fit$convergence != 0) {
    warning(
      "TMB optimization did not converge. ",
      "Estimates may be unreliable. ",
      "Consider checking starting values, model specification, or data quality."
    )
  }
  
  # Reports
  rep  <- TMB::sdreport(obj)
  rep2 <- summary(rep, select = "report")
  
  # Get sensitivities and specificities at the thresholds
  qq <- qnorm(1-(1-conflevel)/2)
  rls1 <- which(rownames(rep2)=="logitSurv1")
  rls0 <- which(rownames(rep2)=="logitSurv0")
  
  if(testdir=="greater"){
    sens <- data.frame(threshold=dat2$threshold,
                       logitSurv1=rep2[rls1,"Estimate"],
                       Std_Error=rep2[rls1,"Std. Error"],
                       CI_Lower=NA,
                       CI_Upper=NA)
    sens$CI_Lower   <- with(sens,logitSurv1-qq*Std_Error)
    sens$CI_Upper   <- with(sens,logitSurv1+qq*Std_Error)
    sens$Sens       <- with(sens,plogis(logitSurv1))
    sens$SensCI_Lower <- with(sens,plogis(CI_Lower))
    sens$SensCI_Upper <- with(sens,plogis(CI_Upper))
  
    spec <- data.frame(threshold=dat2$threshold,
                       logitSurv0=rep2[rls0,"Estimate"],
                       Std_Error=rep2[rls0,"Std. Error"],
                       CI_Lower=NA,
                       CI_Upper=NA)
    spec$CI_Lower   <- with(spec,logitSurv0-qq*Std_Error)
    spec$CI_Upper   <- with(spec,logitSurv0+qq*Std_Error)
    spec$Spec       <- with(spec,1-plogis(logitSurv0))
    spec$SpecCI_Lower <- with(spec,1-plogis(CI_Upper))
    spec$SpecCI_Upper <- with(spec,1-plogis(CI_Lower))
   
    sesp <- data.frame(threshold=sens$threshold,
                       conflevel=conflevel,
                       Sens=sens$Sens,
                       SensCI_Lower=sens$SensCI_Lower,
                       SensCI_Upper=sens$SensCI_Upper,
                       Spec=spec$Spec,
                       SpecCI_Lower=spec$SpecCI_Lower,
                       SpecCI_Upper=spec$SpecCI_Upper)
  }
  ########
  if(testdir=="less"){
    sens <- data.frame(threshold=dat2$threshold,
                       logitSurv1=rep2[rls1,"Estimate"],
                       Std_Error=rep2[rls1,"Std. Error"],
                       CI_Lower=NA,
                       CI_Upper=NA)
    sens$CI_Lower   <- with(sens,logitSurv1-qq*Std_Error)
    sens$CI_Upper   <- with(sens,logitSurv1+qq*Std_Error)
    sens$Sens       <- with(sens,1-plogis(logitSurv1))
    sens$SensCI_Lower <- with(sens,1-plogis(CI_Upper))
    sens$SensCI_Upper <- with(sens,1-plogis(CI_Lower))
    
    spec <- data.frame(threshold=dat2$threshold,
                       logitSurv0=rep2[rls0,"Estimate"],
                       Std_Error=rep2[rls0,"Std. Error"],
                       CI_Lower=NA,
                       CI_Upper=NA)
    spec$CI_Lower   <- with(spec,logitSurv0-qq*Std_Error)
    spec$CI_Upper   <- with(spec,logitSurv0+qq*Std_Error)
    spec$Spec       <- with(spec,plogis(logitSurv0))
    spec$SpecCI_Lower <- with(spec,plogis(CI_Lower))
    spec$SpecCI_Upper <- with(spec,plogis(CI_Upper))
    
    sesp <- data.frame(threshold=sens$threshold,
                       conflevel=conflevel,
                       Sens=sens$Sens,
                       SensCI_Lower=sens$SensCI_Lower,
                       SensCI_Upper=sens$SensCI_Upper,
                       Spec=spec$Spec,
                       SpecCI_Lower=spec$SpecCI_Lower,
                       SpecCI_Upper=spec$SpecCI_Upper)
  }
  # Result object
  res <- list(
    data         = datao,
    restructured = datar,
    fit          = fit,
    sdreport     = rep,
    sdreport2    = rep2,
    distcode     = dat2$dist,
    sensspec     = sesp
  )

  # Assign class
  class(res) <- c("HoyerAFT","Cochrane")

  return(res)
}
