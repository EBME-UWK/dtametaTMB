#' dtametaTMB: Diagnostic Test Accuracy Meta-Analysis using Template Model Builder
#'
#' Functions for fitting diagnostic test accuracy (DTA) meta-analysis models
#' using Template Model Builder (TMB).
#'
#' Implemented methods
#'
#' * Reitsma model
#' * Rutter-Gatsonis HSROC model
#' * Hoyer multiple-threshold model
#' * Latent class analysis for studies with an imperfect reference standard
#'
#' Main workflow
#'
#' * fitReitsma() -> print() -> summary() -> plot() -> forest()
#' * fitReitsmaTMB() -> print() -> summary() -> plot() -> forest()
#' * fitReitsmaSubgroup() -> print() -> summary() -> plot() -> forest()
#' * fitReitsmaSubgroupTMB() -> print() -> summary() -> plot() -> forest()
#' * fitRutterGatsonis() -> print() -> summary() -> plot() -> forest()
#' * fitRutterGatsonisSubgroup() -> print() -> summary() -> plot() -> forest()
#' * fitHoyer() -> print() -> summary() -> plot() -> forest()
#' * fitReitsmaLCA() -> print() -> summary() -> plot() -> forest()
#' * fitReitsmaSubgroupLCA() -> print() -> summary() -> plot() -> forest()
#' * fitRutterGatsonisLCA() -> print() -> summary() -> plot() -> forest()
#' * fitRutterGatsonisSubgroupLCA() -> print() -> summary() -> plot() -> forest()
#'
#' Included datasets
#'
#' * anaemia
#' * anticcp
#' * diabetes
#' * dementia
#' * FENO
#' * pap
#' * RF
#' * schuetz
#' * tub
#'
#' See the package vignettes for worked examples and model descriptions.
#'
#' @name dtametaTMB
#' @rawNamespace useDynLib(dtametaTMB, .registration=TRUE); useDynLib(dtametaTMB_TMBExports)
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL