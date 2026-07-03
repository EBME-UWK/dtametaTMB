#' @export
logLik.RutterGatsonis <- function(object, ...) {
 structure(
    -object$fit$objective,
    class = "logLik",
    df = length(object$fit$par),
    nobs = nrow(object$data)
  )
}

#' @export
logLik.RutterGatsonisSubgroup <- function(object, ...) {
  structure(
    -object$fit$objective,
    class = "logLik",
    df = length(object$fit$par),
    nobs = nrow(object$data)
  )
}

#' @export
logLik.RutterGatsonisReg <- function(object, ...) {
  structure(
    -object$fit$objective,
    class = "logLik",
    df = length(object$fit$par),
    nobs = nrow(object$data)
  )
}

#' @export
logLik.HoyerAFT <- function(object, ...) {
  structure(
    -object$fit$objective,
    class = "logLik",
    df = length(object$fit$par),
    nobs = length(unique(object$data$study))
  )
}

#' @export
logLik.Reitsma <- function(object, ...){
  structure(
    -object$glmmTMB$fit$objective,
    class = "logLik",
    df = length(object$glmmTMB$fit$par),
    nobs = nrow(object$data)
  )
}


#' @export
logLik.ReitsmaSubgroup <- function(object, ...){
  structure(
    -object$glmmTMB_mu$fit$objective,
    class = "logLik",
    df = length(object$glmmTMB_mu$fit$par),
    nobs = nrow(object$data)
  )
}