#' Compute Initial Parameter Values for Threshold-Based bivariate time-to-event DTA Models
#'
#' Estimates starting values for fixed and random effects parameters used
#' in likelihood-based diagnostic test accuracy (DTA) models based on
#' interval-censored data (e.g., Hoyer et al., 2018).
#'
#' The function fits separate intercept-only parametric survival models
#' for diseased and non-diseased groups using weighted interval-censored
#' likelihoods. It further derives initial estimates of between-study
#' variability based on study-specific weighted averages of log-thresholds.
#'
#' @param restructured A data frame in interval format as produced by
#'   \code{\link{restructure}} (specifically the \code{restructured}
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
#' Separate intercept-only survival models are fitted for diseased and
#' non-diseased observations using interval-censored likelihoods weighted
#' by event counts.
#'
#' Initial values for the random effects are obtained from the covariance
#' matrix of study-specific weighted mean log-thresholds for diseased and
#' non-diseased groups.
#'
#' If only a single study is available, a warning is issued as the
#' estimation of between-study variability is not reliable.
#'
#' @examples
#' \dontrun{
#' # Restructure data
#' data("diabetes")
#' res <- restructure(
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
#'
#' # Compute initial parameters
#' init <- getInitParms(res$restructured)
#' }
#'
#' @references
#' Hoyer, A., Hirt, S., Kuss, O. (2018).
#' Meta-analysis of full ROC curves using bivariate time-to-event models
#' for interval-censored data.
#' \emph{Research Synthesis Methods}, 9(1), 1759--2879.
#'
#' @importFrom survival survreg Surv
#' @importFrom stats aggregate cov cov2cor
#' @export
getInitParms <- function(restructured, dist="loglogistic") {
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
  
  datfit0 <- restructured[restructured$events0!=0,]
  datfit1 <- restructured[restructured$events1!=0,]

  fit0 <- survival::survreg(Surv(lowerB, upperB, type = "interval2") ~ 1,
                            data = datfit0,
                            weights = datfit0$events0,
                            dist = dist)

  fit1 <- survival::survreg(Surv(lowerB, upperB, type = "interval2") ~ 1,
                            data = restructured[restructured$events1!=0,],
                            weights = datfit1$events1,
                            dist = dist)

  beta0_init   <- fit0$coefficients
  lambda0_init <- max(fit0$scale,1e-5)
  beta1_init   <- fit1$coefficients
  lambda1_init <- max(fit1$scale,1e-5)

  ## random effects
  lmeantest0 <- aggregate(cbind(lcutmean * restructured$events0, restructured$events0) ~ study,
                          data = restructured,
                          FUN = sum)
  lmeantest0$lmeantest0 <- lmeantest0[, 2] / lmeantest0[, 3]
  lmeantest1 <- aggregate(cbind(lcutmean * restructured$events1, restructured$events1) ~ study,
                          data = restructured,
                          FUN = sum)
  lmeantest1$lmeantest1 <- lmeantest1[, 2] / lmeantest1[, 3]
  meantest <- merge(lmeantest0[c("study", "lmeantest0")],
                    lmeantest1[c("study", "lmeantest1")],
                    by = "study")
  ssc           <- stats::cov(meantest[, c("lmeantest0", "lmeantest1")])
  su0_init      <- max(sqrt(ssc[1,1]),1e-5)
  su1_init      <- max(sqrt(ssc[2,2]),1e-5)
  coru0u1_init  <- min(max(cov2cor(ssc)[1,2],-0.99),0.99)
  if(dist=="weibull"){ distcode = 1 }
  if(dist=="lognormal"){ distcode = 2 }
  if(dist=="loglogistic"){ distcode = 3 }
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
