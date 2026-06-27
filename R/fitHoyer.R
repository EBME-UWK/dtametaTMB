#' Fit Threshold-Based Bivariate Time-to-Event DTA Model (Hoyer AFT)
#'
#' Fits the threshold-based bivariate accelerated failure time (AFT)
#' model for diagnostic test accuracy (DTA) meta-analysis as described
#' by Hoyer et al. (2018). The model is estimated using
#' interval-censored likelihoods via Template Model Builder (TMB).
#'
#' @description
#' This is a high-level wrapper that performs the full workflow for
#' fitting the Hoyer AFT model:
#' \enumerate{
#'   \item Restructuring of threshold-based data into interval format
#'   \item Estimation of initial parameter values
#'   \item Model fitting via likelihood maximization
#' }
#'
#' The input data must contain cumulative counts for true positives,
#' false positives, false negatives, and true negatives at multiple
#' thresholds within each study.
#'
#' @param data A data frame containing study-level diagnostic test data.
#'
#' @param TP,FP,FN,TN Unquoted column names in \code{data} corresponding
#'   to true positives, false positives, false negatives, and true negatives.
#'
#' @param study Unquoted column name identifying studies.
#'
#' @param threshold Unquoted column name specifying the observed threshold values.
#'   Must be positive and strictly increasing within each study.
#'
#' @param smallest Numeric. Positive lower bound defining the leftmost
#'   interval. Must be smaller than the minimum observed threshold.
#'
#' @param largest Numeric. Positive upper bound defining the rightmost
#'   interval. Must be larger than the maximum observed threshold.
#'   
#' @param eval_threshold Optional numeric value or vector specifying the 
#'   prediction grid threshold(s) at which sensitivity and specificity 
#'   should be evaluated. If \code{NA} (default), the median threshold from 
#'   the original data is used.
#'
#' @param dist Character string specifying the parametric distribution
#'   for the AFT model. One of \code{"weibull"}, \code{"lognormal"}, or
#'   \code{"loglogistic"} (default).
#' 
#' @param testdirection Direction of the test. Enter \code{"greater"} 
#'   when larger test values indicate disease. 
#'   Conversely, enter \code{"less"} when lower test values indicate
#'   disease (e.g. anaemia-type tests). 
#'   Defaults to \code{"greater"}.
#'
#' @param conflevel Confidence level for confidence intervals for sensitivities
#'   and specificities at the chosen thresholds. Defaults to \code{0.95}.
#'
#' @param ... Additional arguments passed to \code{\link{fitHoyerAFT}}.
#'
#' @return
#' An object of class \code{"HoyerAFT"} as returned by
#' \code{\link{fitHoyerAFT}}, containing:
#' \describe{
#'   \item{data}{Processed original data}
#'   \item{restructured}{Interval-formatted data}
#'   \item{fit}{Optimization output}
#'   \item{sdreport}{TMB report object}
#'   \item{sdreport2}{Summary of reported parameters}
#'   \item{sensspec}{Sensitivity and specificity estimates}
#' }
#'
#' @details
#' The function internally calls:
#' \itemize{
#'   \item \code{\link{restructure}} to convert cumulative counts into
#'         interval-censored data
#'   \item \code{\link{getInitParms}} to estimate starting values
#'   \item \code{\link{fitHoyerAFT}} to fit the model
#' }
#'
#' Missing values in the input data are removed prior to analysis.
#' A message is issued listing the affected studies.
#'
#' @references
#' Hoyer, A., Hirt, S., Kuss, O. (2018).
#' Meta-analysis of full ROC curves using bivariate time-to-event models for interval-censored data.
#' \emph{Research Synthesis Methods}, 9(1), 62-72.
#' \doi{10.1002/jrsm.1273}
#'
#' @examples
#' \dontrun{
#' data("diabetes")
#'
#' fit <- fitHoyer(
#'   data = diabetes,
#'   TP = TP,
#'   FP = FP,
#'   FN = FN,
#'   TN = TN,
#'   study = study,
#'   threshold = threshold,
#'   testdirection = "greater",
#'   smallest = 2,
#'   largest = 10
#' )
#'
#' summary(fit)
#' }
#'
#' @export

fitHoyer <- function(data,
                     TP, FP, FN, TN,
                     study,
                     threshold,
                     smallest,
                     largest,
                     testdirection=c("greater","less"),
                     eval_threshold=NA,
                     conflevel=0.95,
                     dist = "loglogistic",
                     ...) {
  
  # Basic argument checks
  if (missing(threshold)) {
    stop("'threshold' must be provided.")
  }
  if (missing(smallest) || missing(largest)) {
    stop("'smallest' and 'largest' must be provided.")
  }
  
  # Step 1: Restructure data
  res <- restructure(
    data = data,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    study = study,
    threshold = threshold,
    testdirection = testdirection,
    smallest = smallest,
    largest = largest
  )
  
  # Step 2: Initial parameters
  init <- getInitParms(res$restructured, dist = dist)
  
  # Step 3: Fit model
  fit <- fitHoyerAFT(res, init, threshold=eval_threshold, conflevel=conflevel)
  
  return(fit)
}
