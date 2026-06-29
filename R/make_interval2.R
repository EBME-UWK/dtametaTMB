#' Prepare Interval-Censored Data for the Hoyer AFT Model
#'
#' Internal helper function that restructures diagnostic test accuracy (DTA)
#' data with multiple thresholds per study into the interval-censored format
#' required by the Hoyer accelerated failure time (AFT) model.
#'
#' @details
#' For studies with multiple thresholds, the function constructs consecutive
#' censoring intervals between ordered threshold values. The first interval
#' represents left censoring, intermediate intervals represent interval
#' censoring between adjacent thresholds, and the final interval represents
#' right censoring.
#'
#' For each interval, cumulative counts are transformed into interval-specific
#' event counts for diseased and non-diseased subjects. These counts form the
#' basis of the interval-censored likelihood used by the Hoyer AFT model.
#'
#' The arguments `smallest` and `largest` define finite lower and upper
#' reference values used to approximate the mean log-threshold (`lcutmean`)
#' for the boundary intervals.
#'
#' Studies with a single threshold are handled automatically by constructing
#' the corresponding boundary intervals.
#'
#' @param dat A data frame containing the restructured data for a single study.
#' @param smallest Positive value used as the lower boundary for the
#'   left-censored interval.
#' @param largest Value used as the upper boundary for the right-censored
#'   interval.
#'
#' @return
#' A data frame containing one row per censoring interval with the interval
#' boundaries, interval-specific event counts, censoring type, and auxiliary
#' quantities required for likelihood evaluation.
#'
#' @seealso
#' \code{\link{restructure_data}},
#' \code{\link{initHoyer}},
#' \code{\link{fitHoyer}},
#' \code{\link{fitHoyerAFT}}
#'
#' @keywords internal
#' @noRd
make_interval2 <- function(dat,smallest,largest) {
  dat <- dat[order(dat$threshold), ]
  n   <- nrow(dat)
  out <- vector("list", n + 1)
  ## First interval
  out[[1]] <- data.frame(
    study = dat$study[1],
    TP = dat$TP[1],
    TN = dat$TN[1],
    D = dat$D[1],
    H = dat$H[1],
    threshold = dat$threshold[1],
    lowerB = NA,
    upperB = dat$threshold[1],
    events1 = dat$FN[1],
    events0 = dat$TN[1],
    ctype = 1,
    lcutmean = (log(dat$threshold[1])+log(smallest))/2
  )
  ## Intermediate intervals
  if (n > 1) {
    for (i in 2:n) {
      out[[i]] <- data.frame(
        study = dat$study[i],
        TP = dat$TP[i],
        TN = dat$TN[i],
        D = dat$D[i],
        H = dat$H[i],
        threshold = dat$threshold[i],
        lowerB = dat$threshold[i - 1],
        upperB = dat$threshold[i],
        events1 = dat$FN[i] - dat$FN[i - 1],
        events0 = dat$TN[i] - dat$TN[i - 1],
        ctype = 2,
        lcutmean = (log(dat$threshold[i-1])+log(dat$threshold[i]))/2
      )
    }
  }
  ## Final right-censored interval
  out[[n + 1]] <- data.frame(
    study = dat$study[n],
    TP = dat$TP[n],
    TN = dat$TN[n],
    D = dat$D[n],
    H = dat$H[n],
    threshold = dat$threshold[n],
    lowerB = dat$threshold[n],
    upperB = NA,
    events1 = dat$TP[n],
    events0 = dat$FP[n],
    ctype = 3,
    lcutmean = (log(largest)+log(dat$threshold[n]))/2
  )
  do.call(rbind, out)
}
