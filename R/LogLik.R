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


#' @importFrom stats logLik pchisq
NULL

anova_lrt <- function(object, ..., test = "Chisq") {
  models <- list(object, ...)
  if (length(models) < 2) {
    stop("At least two fitted models are required.")
  }
  ll          <- lapply(models, stats::logLik)
  logLik_vals <- sapply(ll, as.numeric)
  dfs         <- sapply(ll, attr, which = "df")
  ord         <- order(dfs)
  logLik_vals <- logLik_vals[ord]
  dfs         <- dfs[ord]
  if(any(diff(dfs) == 0)) {
    warning(
      "Likelihood-ratio tests require nested models with different numbers of parameters.",
      "Likelihood-ratio tests may not be meaningful."
    )
  }
  models      <- models[ord]
  Df          <- c(NA, diff(dfs))
  Chisq       <- c(NA, 2 * diff(logLik_vals))
  if(any(Chisq[-1] < -1e-8)) {
    warning(
      "Likelihood decreased when moving to a more complex model. ",
      "Models may not be nested."
    )
  }
  `Pr(>Chisq)`<- c(NA,stats::pchisq(Chisq[-1], df = Df[-1], lower.tail = FALSE))
  out <- data.frame(Df      = dfs,
                    logLik  = logLik_vals,
                    Df.diff = Df,
                      Chisq = Chisq,
               `Pr(>Chisq)` = `Pr(>Chisq)`,
                check.names = FALSE)
  
  rownames(out) <- paste0("Model ", seq_len(nrow(out)))
  class(out)    <- c("anova", "data.frame")
  return(out)
}

#' @export
anova.Reitsma <- function(object, ..., test = "Chisq") {
  anova_lrt(object, ..., test = test)
}

#' @export
anova.ReitsmaSubgroup <- function(object, ..., test = "Chisq") {
  anova_lrt(object, ..., test = test)
}

#' @export
anova.RutterGatsonis <- function(object, ..., test = "Chisq") {
  anova_lrt(object, ..., test = test)
}

#' @export
anova.RutterGatsonisSubgroup <- function(object, ..., test = "Chisq") {
  anova_lrt(object, ..., test = test)
}

#' @export
anova.RutterGatsonisReg <- function(object, ..., test = "Chisq") {
  anova_lrt(object, ..., test = test)
}



