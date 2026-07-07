#' Construct Interval Data for Threshold-Based Bivariate Time-to-Event Models
#'
#' Transforms study-level diagnostic test accuracy (DTA) data with multiple thresholds
#' into interval-censored format suitable for likelihood-based modelling
#' (Hoyer et al., 2018).
#'
#' The function constructs:
#' \itemize{
#'   \item A left-censored interval (below the first threshold)
#'   \item Intermediate intervals between adjacent thresholds
#'   \item A right-censored interval (above the last threshold)
#' }
#'
#' Event counts are derived from cumulative counts to represent the number
#' of observations falling within each interval for diseased and
#' non-diseased groups. Therefore, input counts must be cumulative over
#' increasing thresholds within each study.
#'
#' @param data A data.frame containing study-level data.
#' @param TP True positives (column name).
#' @param FP False positives (column name).
#' @param FN False negatives (column name).
#' @param TN True negatives (column name).
#' @param study Study identifier (column name).
#' @param threshold The observed threshold value (column name).
#'   Must be positive and strictly increasing within each study.
#'
#' @param testdirection Direction of the test. Enter \code{"greater"} 
#'   when larger test values indicate disease. 
#'   Conversely, enter \code{"less"} when lower test values indicate
#'   disease (e.g. anaemia-type tests). 
#'   Defaults to \code{"greater"}.
#'
#' @param smallest Positive lower bound used to define the
#'   leftmost interval. Must be smaller than minimum \code{threshold}.
#'
#' @param largest Positive upper bound used to define the
#'   rightmost interval. Must be greater than maximum \code{threshold}.
#' 
#'
#' @return A list with two components:
#' \describe{
#'   \item{restructured}{A data frame with one row per constructed interval containing:
#'     \describe{
#'       \item{study}{Study identifier}
#'       \item{TP, TN, n1, n0}{Original counts carried forward}
#'       \item{threshold}{Threshold associated with the interval}
#'       \item{lowerB}{Lower interval bound (NA for left-censored)}
#'       \item{upperB}{Upper interval bound (NA for right-censored)}
#'       \item{events1}{Number of diseased observations in the interval}
#'       \item{events0}{Number of non-diseased observations in the interval}
#'       \item{ctype}{Censoring type (1 = left, 2 = interval, 3 = right)}
#'       \item{lcutmean}{Midpoint of log-thresholds defining the interval}
#'     }
#'   }
#'   \item{original}{The processed original data including (derived) quantities:
#'     \describe{
#'       \item{n1}{Total number of diseased individuals (TP + FN)}
#'       \item{n0}{Total number of non-diseased individuals (TN + FP)}
#'       \item{sens}{Sensitivity (TP / n1)}
#'       \item{spec}{Specificity (TN / n0)}
#'       \item{fpr}{False positive rate (FP / n0)}
#'       \item{testdirection}{As specified by \code{testdirection}}
#'     }
#'   }
#' }
#'
#' @details
#' The function first validates that counts are numeric,
#' non-negative integers and that threshold values are positive.
#'
#' Within each study, total numbers of diseased (\code{n1}) and
#' non-diseased (\code{n0}) individuals are required to be constant across
#' thresholds. Intermediate interval counts are computed as differences
#' between cumulative counts.
#'
#' The log-scale midpoint (\code{lcutmean}) is provided to support
#' initialization of random effects in subsequent Hoyer AFT models.
#'
#' For \code{testdirection = "greater"}, the function assumes that sensitivity 
#' decreases and specificity increases with increasing thresholds. For 
#' \code{testdirection = "less"}, the reverse monotonicity is required.
#' Internally, the function standardizes the definition of a positive test 
#' result so that it always corresponds to values above the threshold. For 
#' \code{testdirection = "less"}, this is achieved by relabeling the observed 
#' counts. This does not affect the resulting sensitivity, specificity, or ROC 
#' curve, but ensures compatibility with the model formulation.
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
#'
#' @references
#' Hoyer, A., Hirt, S., Kuss, O. (2018).
#' Meta-analysis of full ROC curves using bivariate time-to-event models for interval-censored data.
#' \emph{Research Synthesis Methods}, 9(1), 62-72.
#' \doi{10.1002/jrsm.1273}
#' 
#' @importFrom stats complete.cases
#' @export
restructure_data <- function(data,
                             TP, FP, FN, TN,
                             threshold,
                             study,
                             smallest,
                             largest,
                             testdirection=c("greater", "less")) {
  
    testdirection <- match.arg(testdirection)
  
    if (!is.data.frame(data)) {
      stop("'data' must be a data.frame.")
    }
    
    TP_col <- deparse(substitute(TP))
    FP_col <- deparse(substitute(FP))
    FN_col <- deparse(substitute(FN))
    TN_col <- deparse(substitute(TN))
    threshold_col <- deparse(substitute(threshold))
    study_col <- deparse(substitute(study))
    
    dat <- data.frame(
      study = data[[study_col]],
      TP = data[[TP_col]],
      TN = data[[TN_col]],
      FP = data[[FP_col]],
      FN = data[[FN_col]],
      threshold = data[[threshold_col]]
    )
    
  excluded <- !stats::complete.cases(dat)
  if (any(excluded)) {
    removed_studies <- unique(dat$study[excluded])
    message(
      "Removed rows with missing values for studies: ",
      paste(removed_studies, collapse = ", ")
    )
  }
  
  dat <- dat[stats::complete.cases(dat),]

  # Validation
  numeric_cols <- c("TP", "TN", "FP", "FN", "threshold")
  non_numeric <- numeric_cols[!sapply(dat[numeric_cols], is.numeric)]
  if (length(non_numeric) > 0) {
    stop("Columns must be numeric: ", paste(non_numeric, collapse = ", "))
  }

  # Columns
  count_cols <- c("TP", "TN", "FP", "FN")
  # Check for non-integers or negative values
  invalid_counts <- sapply(dat[count_cols], function(x) {
    any(x < 0 | x != floor(x), na.rm = TRUE)
  })

  if (any(invalid_counts)) {
    stop("Columns TP, TN, FP, FN must contain non-negative integer counts.")
  }

  if (any(dat$threshold <= 0, na.rm = TRUE)) {
    stop("'threshold' must be positive.")
  }

  if (!is.numeric(smallest) || smallest <= 0) {
    stop("'smallest' must be positive.")
  }

  if (!is.numeric(largest) || largest <= 0) {
    stop("'largest' must be positive.")
  }

  if (smallest >= min(dat$threshold,na.rm=TRUE)){
    stop("'smallest' must be smaller than minimum threshold.")
  }
  
  if (largest <= max(dat$threshold,na.rm=TRUE)){
    stop("'largest' must be larger than maximum threshold.")
  }
  
  # Order according to study and reported threshold
  dat <- dat[order(dat$study, dat$threshold), ]
  # Derived totals
  dat$n1 <- dat$TP + dat$FN
  dat$n0 <- dat$TN + dat$FP

  check_consistency <- function(df) {
    if (length(unique(df$n1)) > 1) {
      stop(
        sprintf(
          paste0(
            "Inconsistent diseased counts (TP+FN) within study '%s'.\n",
            "Rows: %s\n",
            "Counts: %s"
          ),
          df$study[1],
          paste(rownames(df), collapse = ", "),
          paste(unique(df$n1), collapse = ", ")
        )
      )
    }
    if (length(unique(df$n0)) > 1) {
      stop(
        sprintf(
          paste0(
            "Inconsistent non-diseased counts (TN+FP) within study '%s'.\n",
            "Rows: %s\n",
            "Counts: %s"
          ),
          df$study[1],
          paste(rownames(df), collapse = ", "),
          paste(unique(df$n0), collapse = ", ")
        )
      )
    }
    TRUE
  }
  # Shared helper: build a properly separated, multi-line message
  # for one or more monotonicity violations detected via which(diff(...) ...).
  .mono_violation_msg <- function(df, col, idx, direction) {
    msgs <- sprintf(
      paste0(
        "%s must be %s with increasing thresholds in study '%s'.\n",
        "  Rows: %s -> %s\n",
        "  Thresholds: %s -> %s\n",
        "  %s: %s -> %s"
      ),
      col, direction, df$study[1],
      rownames(df)[idx], rownames(df)[idx + 1],
      df$threshold[idx], df$threshold[idx + 1],
      col, df[[col]][idx], df[[col]][idx + 1]
    )
    paste(msgs, collapse = "\n\n")
  }
  
  # Check consistency within each study
  check_monotonicity_greater <- function(df) {
    df <- df[order(df$threshold), ]
    # FN non-decreasing
    idx <- which(diff(df$FN) < 0)
    if (length(idx) > 0) {
      stop(.mono_violation_msg(df, "FN", idx, "non-decreasing"), call. = FALSE)
    }
    # TN non-decreasing
    idx <- which(diff(df$TN) < 0)
    if (length(idx) > 0) {
      stop(.mono_violation_msg(df, "TN", idx, "non-decreasing"), call. = FALSE)
    }
    # TP non-increasing
    idx <- which(diff(df$TP) > 0)
    if (length(idx) > 0) {
      stop(.mono_violation_msg(df, "TP", idx, "non-increasing"), call. = FALSE)
    }
    # FP non-increasing
    idx <- which(diff(df$FP) > 0)
    if (length(idx) > 0) {
      stop(.mono_violation_msg(df, "FP", idx, "non-increasing"), call. = FALSE)
    }
    TRUE
  }
  
  check_monotonicity_less <- function(df) {
    df <- df[order(df$threshold), ]
    # FN non-increasing
    idx <- which(diff(df$FN) > 0)
    if (length(idx) > 0) {
      stop(.mono_violation_msg(df, "FN", idx, "non-increasing"), call. = FALSE)
    }
    # TN non-increasing
    idx <- which(diff(df$TN) > 0)
    if (length(idx) > 0) {
      stop(.mono_violation_msg(df, "TN", idx, "non-increasing"), call. = FALSE)
    }
    # TP non-decreasing
    idx <- which(diff(df$TP) < 0)
    if (length(idx) > 0) {
      stop(.mono_violation_msg(df, "TP", idx, "non-decreasing"), call. = FALSE)
    }
    # FP non-decreasing
    idx <- which(diff(df$FP) < 0)
    if (length(idx) > 0) {
      stop(.mono_violation_msg(df, "FP", idx, "non-decreasing"), call. = FALSE)
    }
    TRUE
  }
  
  # Apply checks per study
  invisible(lapply(split(dat, dat$study), check_consistency))
  
  if(testdirection=="greater"){
    invisible(lapply(split(dat, dat$study), check_monotonicity_greater))
  }
  
  if(testdirection=="less"){
    invisible(lapply(split(dat, dat$study), check_monotonicity_less))
  }
  
  dat$testdirection <- testdirection
  dat$sens <- dat$TP/dat$n1
  dat$spec <- dat$TN/dat$n0
  dat$fpr  <- dat$FP/dat$n0
  dat2     <- dat  

  if (testdirection == "less") {
    tmp_TP <- dat2$TP
    tmp_FP <- dat2$FP
    
    dat2$TP <- dat2$FN
    dat2$FP <- dat2$TN
    dat2$FN <- tmp_TP
    dat2$TN <- tmp_FP
  }
  
  # Apply interval construction
  result <- do.call(
    rbind,
    lapply(split(dat2, dat2$study),
           function(x) make_interval2(x, smallest, largest))
  )
  
  return(list(restructured=result,original=dat))
}
