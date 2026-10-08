#' Compute Initial Parameter Values for HoyerAFT Models
#'
#' 
#' Computes initial values for the threshold-based bivariate time-to-event model
#' of Hoyer et al. (2018).
#'
#' The function fits separate intercept-only parametric survival models
#' for diseased and non-diseased groups using weighted interval-censored
#' likelihoods. It further derives initial estimates of between-study
#' variability based on study-specific weighted averages of log-thresholds.
#'
#' @param restructured A data frame in interval format as produced by
#'   \code{\link{restructure_data}} (specifically the \code{restructured}
#'   component of its output). Must contain the columns:
#'   \describe{
#'     \item{study}{Study identifier}
#'     \item{lowerB}{Lower interval bound}
#'     \item{upperB}{Upper interval bound}
#'     \item{events0}{Non-diseased counts within interval}
#'     \item{events1}{Diseased counts within interval}
#'     \item{ctype}{Censoring type (1 = left, 2 = interval, 3 = right)}
#'     \item{lcutmean}{Midpoint of log-threshold interval}
#'   }
#'
#' @param dist Character string specifying the parametric distribution
#'   used in the survival regression models. Must be one of:
#'   \code{"weibull"}, \code{"lognormal"}, or \code{"loglogistic"}.
#'   Default is \code{"loglogistic"}.
#'
#' @return A single-row data frame containing initial parameter values:
#' \describe{
#'   \item{beta0_init}{Intercept for non-diseased group}
#'   \item{lambda0_init}{Scale parameter for non-diseased group}
#'   \item{beta1_init}{Intercept for diseased group}
#'   \item{lambda1_init}{Scale parameter for diseased group}
#'   \item{su0_init}{Standard deviation of random effects (non-diseased)}
#'   \item{su1_init}{Standard deviation of random effects (diseased)}
#'   \item{coru0u1_init}{Correlation between random effects}
#'   \item{distcode}{Numeric code for the distribution
#'     (1 = Weibull, 2 = lognormal, 3 = loglogistic)}
#' }
#'
#' @details
#' To ensure compatibility with \code{survival::survreg}, interval bounds
#' are modified as follows:
#' \itemize{
#'   \item Left-censored intervals (\code{ctype = 1}) are assigned a small
#'         positive lower bound.
#'   \item Right-censored intervals (\code{ctype = 3}) are assigned an
#'         infinite upper bound.
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
#' 
#'
#' @references
#' Hoyer, A., Hirt, S., Kuss, O. (2018).
#' Meta-analysis of full ROC curves using bivariate time-to-event models
#' for interval-censored data.
#' \emph{Research Synthesis Methods}, 9(1), 62-72.
#' \doi{10.1002/jrsm.1273}
#'
#' @importFrom survival survreg Surv
#' @importFrom stats aggregate sd cor weighted.mean
#' @export
initHoyerAFT <- function(restructured, dist="loglogistic") {
  # Check distribution
  valid_dists <- c("weibull", "lognormal", "loglogistic")
  if (!dist %in% valid_dists) {
    stop("Argument 'dist' must be one of: 'weibull', 'lognormal', 'loglogistic'.")
  }

  # Check required columns
  required <- c("study", "lowerB", "upperB", "events0", "events1", "ctype", "lcutmean")
  missing <- setdiff(required, names(restructured))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }

  # Check non-zero data for model fitting
  if (all(restructured$events0 == 0)) {
    stop("No non-diseased events (events0) available.")
  }
  if (all(restructured$events1 == 0)) {
    stop("No diseased events (events1) available.")
  }

  # Optional: study count
  if (length(unique(restructured$study)) < 2) {
    warning("Only one study: random effects may be unstable.")
  }

  restructured$lowerB[restructured$ctype == 1] <- 1e-09
  restructured$upperB[restructured$ctype == 3] <- Inf
  
  datfit0 <- restructured[restructured$events0>0,,drop=FALSE]
  datfit1 <- restructured[restructured$events1>0,,drop=FALSE]

  fit0 <- tryCatch(suppressWarnings(survival::survreg(survival::Surv(lowerB, upperB, type = "interval2") ~ 1,
                                                      data = datfit0,
                                                      weights = datfit0$events0,
                                                      dist = dist)), error = function(e) NULL)

  fit1 <- tryCatch(suppressWarnings(survival::survreg(survival::Surv(lowerB, upperB, type = "interval2") ~ 1,
                                                      data = datfit1,
                                                      weights = datfit1$events1,
                                                      dist = dist)), error = function(e) NULL)

  beta0_init   <- if(is.null(fit0)) { NA_real_ } else { unname(fit0$coefficients[1]) }
  if (!is.finite(beta0_init)) { beta0_init <- stats::weighted.mean(restructured$lcutmean,w=restructured$events0,na.rm=TRUE) }
  if (!is.finite(beta0_init)) { beta0_init <- 0 }
  ##
  lambda0_init <- if (is.null(fit0)) { NA_real_ } else { fit0$scale }
  if (!is.finite(lambda0_init) || lambda0_init <= 0) { 
    m0 <- stats::weighted.mean(restructured$lcutmean,w=restructured$events0,na.rm=TRUE)
    x0 <- (restructured$lcutmean-m0)**2
    lambda0_init <- sqrt(stats::weighted.mean(x0,w=restructured$events0,na.rm=TRUE)) }
  if (!is.finite(lambda0_init) || lambda0_init <= 0) { lambda0_init <- 1 }
  lambda0_init <- max(lambda0_init, 1e-5)
  ##
  beta1_init   <- if(is.null(fit1)) { NA_real_ } else { unname(fit1$coefficients[1]) }
  if (!is.finite(beta1_init)) { beta1_init <- stats::weighted.mean(restructured$lcutmean,w=restructured$events1,na.rm=TRUE) }
  if (!is.finite(beta1_init)) { beta1_init <- 0 }
  ##
  lambda1_init <- if (is.null(fit1)) { NA_real_ } else { fit1$scale }
  if (!is.finite(lambda1_init) || lambda1_init <= 0) { 
    m1 <- stats::weighted.mean(restructured$lcutmean,w=restructured$events1,na.rm=TRUE)
    x1 <- (restructured$lcutmean-m1)**2
    lambda1_init <- sqrt(stats::weighted.mean(x1,w=restructured$events1,na.rm=TRUE)) }
  if (!is.finite(lambda1_init) || lambda1_init <= 0) { lambda1_init <- 1 }
  lambda1_init <- max(lambda1_init, 1e-5)
  ## random effects
  lmeantest0 <- aggregate(cbind(lcutmean * restructured$events0, restructured$events0) ~ study,
                          data = restructured,
                          FUN = sum)
  lmeantest0$lmeantest0 <- lmeantest0[, 2] / lmeantest0[, 3]
  lmeantest0 <- lmeantest0[is.finite(lmeantest0$lmeantest0),c("study", "lmeantest0"),drop = FALSE]
  #
  lmeantest1 <- aggregate(cbind(lcutmean * restructured$events1, restructured$events1) ~ study,
                          data = restructured,
                          FUN = sum)
  lmeantest1$lmeantest1 <- lmeantest1[, 2] / lmeantest1[, 3]
  lmeantest1 <- lmeantest1[is.finite(lmeantest1$lmeantest1),c("study", "lmeantest1"),drop = FALSE]
  #
  meantest <- merge(lmeantest0,lmeantest1,by="study", all=FALSE,sort=FALSE)
  # Neutral fallbacks
  su0_init <- 0.5
  su1_init <- 0.5
  coru0u1_init <- 0
  #
  if (nrow(lmeantest0) >= 2L) { 
    sd0 <- stats::sd(lmeantest0$lmeantest0)
    if (is.finite(sd0)) {su0_init <- max(sd0, 1e-5)}
  }
  if (nrow(lmeantest1) >= 2L) { 
    sd1 <- stats::sd(lmeantest1$lmeantest1)
    if (is.finite(sd1)) {su1_init <- max(sd1, 1e-5)}
  }
  if (nrow(meantest) >= 3L) {
    sd0_cor <- stats::sd(meantest$lmeantest0)
    sd1_cor <- stats::sd(meantest$lmeantest1)
    if (is.finite(sd0_cor) && is.finite(sd1_cor) && sd0_cor > sqrt(.Machine$double.eps) && sd1_cor > sqrt(.Machine$double.eps)) {
      empirical_cor <- suppressWarnings(stats::cor(meantest$lmeantest0,meantest$lmeantest1))
      if (is.finite(empirical_cor)) {coru0u1_init <- pmin(pmax(empirical_cor, -0.99),0.99)}
    }
  }
  ###
  distcode <- switch(dist,weibull = 1,lognormal = 2,loglogistic = 3)
  res           <- data.frame(beta0_init,
                              lambda0_init,
                              beta1_init,
                              lambda1_init,
                              su0_init,
                              su1_init,
                              coru0u1_init,
                              distcode)
  rownames(res) <- "Initial value"
  return(res)
}
