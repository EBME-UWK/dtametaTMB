#' @keywords internal
#' @noRd
check_preprocess_data <- function(data,TP,FP,FN,TN,study,conflevel){
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame.")
  }
  TP_col <- deparse(substitute(TP))
  FP_col <- deparse(substitute(FP))
  FN_col <- deparse(substitute(FN))
  TN_col <- deparse(substitute(TN))
  study_col <- deparse(substitute(study))

  dat <- data.frame(
    study = data[[study_col]],
    TP = data[[TP_col]],
    TN = data[[TN_col]],
    FP = data[[FP_col]],
    FN = data[[FN_col]])

  excluded <- !stats::complete.cases(dat)
  if (any(excluded)) {
    removed_studies <- unique(dat$study[excluded])
    message(
      "Removed rows with missing values for studies: ",
      paste(removed_studies, collapse = ", ")
    )
  }

  dat <- dat[stats::complete.cases(dat), ]
  # Validation
  numeric_cols <- c("TP", "TN", "FP", "FN")
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
  
  if (!is.numeric(conflevel) || length(conflevel) != 1L ||
      conflevel <= 0 || conflevel >= 1) {
    stop("conflevel must be a single number in (0, 1).")
  }
  
  return(dat)
}