#' HPV Dataset
#'
#' This is the HPV dataset as given in the supplementary material from Nyaga and Arbyn (2023) 
#' from a network meta-analysis of HPV testing and cervical cytology for the detection of CIN2+.
#'
#' @format A data frame with 36 rows and 9 variables:
#' \describe{
#'   \item{StudyID}{Study identifier}
#'   \item{Year}{Year of publication or study}
#'   \item{Age}{Age range of the study population}
#'   \item{Location}{Geogarphical location of the study}
#'   \item{Test}{Which test was used, Hybrid capture 2 (HC2), liquid-based cytology (LBC), or conventional cytology (CC)}
#'   \item{TP}{Number of true positives}
#'   \item{FP}{Number of false positives}
#'   \item{TN}{Number of true negatives}
#'   \item{FN}{Number of false negatives}
#'
#' }
#' @source Nyaga, K., & Arbyn, M. (2023).
#' *Comparison and validation of metadta for meta-analysis of diagnostic test accuracy studies.*
#' Research Synthesis Methods, 14(3), 544-562.
#' \doi{10.1002/jrsm.1634}
"HPV"
