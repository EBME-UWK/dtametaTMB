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
#' * fitReitsma() -> plot() -> forest()
#' * fitReitsmaSubgroup() -> plot() -> forest()
#' * fitRutterGatsonis() -> plot() -> forest()
#' * fitRutterGatsonisSubgroup() -> plot() -> forest()
#' * fitHoyer() -> plot() -> forest()
#' * fitReitsmaLCA() -> plot() -> forest()
#' * fitReitsmaSubgroupLCA() -> plot() -> forest()
#' * fitRutterGatsonisLCA() -> plot() -> forest()
#' * fitRutterGatsonisSubgroupLCA() -> plot() -> forest()
#'
#' Included datasets
#'
#' * anaemia
#' * anticcp
#' * diabetes
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