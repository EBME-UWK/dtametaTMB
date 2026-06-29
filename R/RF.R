#' RF dataset
#'
#' This is the RF data set as given by Nishimura (2007) from a review of rheumatoid factor (RF) to diagnose rheumatoid arthritis.
#'
#' @format A data frame with 37 rows and 7 variables:
#' \describe{
#'   \item{study}{Study identifier}
#'   \item{year}{Year of the study}
#'   \item{TP}{Number of true positives}
#'   \item{FP}{Number of false positives}
#'   \item{FN}{Number of false negatives}
#'   \item{TN}{Number of true negatives}
#'   \item{cutoff}{Threshold used in U/mL for test positivity}
#'   \item{method}{Different methods used to perform the test}
#'
#' }
#' @source Nishimura, K. et al. (2007).
#' *Meta-analysis: diagnostic accuracy of anti-cyclic citrullinated peptide antibody and rheumatoid factor for rheumatoid arthritis*.
#' Annals of Internal Medicine, 146(11), 392-403.
#' \doi{10.7326/0003-4819-146-11-200706050-00008}
"RF"
