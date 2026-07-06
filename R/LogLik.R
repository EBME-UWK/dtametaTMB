#' Extract Log-Likelihood from Diagnostic Test Accuracy Models
#'
#' Returns the maximized log-likelihood for a fitted diagnostic test
#' accuracy meta-analysis model. The returned value is an object of class
#' `"logLik"` and includes the number of estimated parameters (`df`) and
#' the number of observations used in model fitting (`nobs`).
#'
#' @param object A fitted model object.
#' @param ... Not used.
#'
#' @return An object of class `"logLik"`.
#'
#' @seealso
#' [anova.dtametaTMB()]
#'
#' @name logLik.dtametaTMB


#' @rdname logLik.dtametaTMB
#' @export
logLik.RutterGatsonis <- function(object, ...) {
 structure(
    -object$fit$objective,
    class = "logLik",
    df = length(object$fit$par),
    nobs = nrow(object$data)
  )
}

#' @rdname logLik.dtametaTMB
#' @export
logLik.RutterGatsonisSubgroup <- function(object, ...) {
  structure(
    -object$fit$objective,
    class = "logLik",
    df = length(object$fit$par),
    nobs = nrow(object$data)
  )
}

#' @rdname logLik.dtametaTMB
#' @export
logLik.RutterGatsonisReg <- function(object, ...) {
  structure(
    -object$fit$objective,
    class = "logLik",
    df = length(object$fit$par),
    nobs = nrow(object$data)
  )
}

#' @rdname logLik.dtametaTMB
#' @export
logLik.HoyerAFT <- function(object, ...) {
  structure(
    -object$fit$objective,
    class = "logLik",
    df = length(object$fit$par),
    nobs = length(unique(object$data$study))
  )
}

#' @rdname logLik.dtametaTMB
#' @export
logLik.Reitsma <- function(object, ...){
  structure(
    -object$glmmTMB$fit$objective,
    class = "logLik",
    df = length(object$glmmTMB$fit$par),
    nobs = nrow(object$data)
  )
}

#' @rdname logLik.dtametaTMB
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

#' Likelihood Ratio Tests for Diagnostic Test Accuracy Models
#'
#' Compare nested diagnostic test accuracy meta-analysis models using
#' likelihood ratio tests (LRTs).
#'
#' These methods compute likelihood ratio statistics from the fitted model
#' log-likelihoods and compare models ordered by increasing numbers of
#' estimated parameters. The test statistic is
#'
#' \deqn{2(\ell_1 - \ell_0)}
#'
#' where \eqn{\ell_1} and \eqn{\ell_0} are the log-likelihoods of two nested
#' models. Under standard regularity conditions, the statistic follows a
#' chi-squared distribution with degrees of freedom equal to the difference
#' in the numbers of estimated parameters.
#'
#' The models supplied should be nested and fitted to the same dataset.
#' Warnings are issued when models have identical numbers of parameters or
#' when the log-likelihood decreases for a more complex model, suggesting
#' that the nesting assumptions may be violated.
#'
#' @param object A fitted model object.
#' @param ... Additional fitted model objects to be compared.
#' @param test Character string specifying the test to perform. Currently,
#'   only `"Chisq"` is supported.
#'
#' @return An object of class `"anova"` inheriting from `"data.frame"` with
#'   the following columns:
#'   \describe{
#'     \item{Df}{Number of estimated parameters in the model.}
#'     \item{logLik}{Model log-likelihood.}
#'     \item{Df.diff}{Difference in parameters compared with the previous
#'       model.}
#'     \item{Chisq}{Likelihood ratio chi-squared statistic.}
#'     \item{Pr(>Chisq)}{P-value from the chi-squared test.}
#'   }
#'
#' @seealso
#' [logLik()], [anova()]
#'
#' @name anova.dtametaTMB


#' @rdname anova.dtametaTMB
#' @export
anova.Reitsma <- function(object, ..., test = "Chisq") {
  anova_lrt(object, ..., test = test)
}

#' @rdname anova.dtametaTMB
#' @export
anova.ReitsmaSubgroup <- function(object, ..., test = "Chisq") {
  anova_lrt(object, ..., test = test)
}

#' @rdname anova.dtametaTMB
#' @export
anova.RutterGatsonis <- function(object, ..., test = "Chisq") {
  anova_lrt(object, ..., test = test)
}

#' @rdname anova.dtametaTMB
#' @export
anova.RutterGatsonisSubgroup <- function(object, ..., test = "Chisq") {
  anova_lrt(object, ..., test = test)
}

#' @rdname anova.dtametaTMB
#' @export
anova.RutterGatsonisReg <- function(object, ..., test = "Chisq") {
  anova_lrt(object, ..., test = test)
}


