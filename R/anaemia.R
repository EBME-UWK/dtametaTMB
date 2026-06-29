#' Anaemia (synthetic data set)
#' 
#' This is a synthetic data set where haemoglobin was measured 
#' by a point-of-care device, with laboratory measurement acting as 
#' the reference standard. Lower values indicate disease (anaemia).
#'
#' @format A data frame with 81 rows and 11 variables:
#' \describe{
#'   \item{study}{Study identifier}
#'   \item{threshold}{Haemoglobin in g/dL}
#'   \item{TP}{Number of true positives}
#'   \item{FN}{Number of false negatives}
#'   \item{FP}{Number of false positives}
#'   \item{TN}{Number of true negatives}
#'   \item{D}{Number of diseased individuals}
#'   \item{H}{Number of healthy individuals}
#'   \item{sens}{Sensitivity}
#'   \item{fpr}{False positive rate}
#'   \item{testdirection}{Smaller values indicate disease (=less)}
#' }
#' @source Synthetic dataset.
"anaemia"