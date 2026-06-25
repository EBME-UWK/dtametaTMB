#' Print Hoyer AFT Model Object
#'
#' Displays a concise summary of a fitted Hoyer AFT model, including
#' distribution, number of studies, threshold, and optimization status.
#'
#' @param x An object of class \code{"HoyerAFT"}.
#' @param ... Further arguments (unused).
#'
#' @return The input object, returned invisibly.
#'
#' @method print HoyerAFT
#' @export
print.HoyerAFT <- function(x, ...) {

  cat("\nHoyer AFT Model Fit\n")
  cat("-------------------\n")

  # Distribution
  dist_name <- switch(as.character(x$distcode),
                      "1" = "Weibull",
                      "2" = "Lognormal",
                      "3" = "Loglogistic",
                      "Unknown")

  cat("Distribution :", dist_name, "\n")

  # Number of studies
  nstudy <- length(unique(x$restructured$study))
  cat("Number of studies :", nstudy, "\n")
  
  # Test direction
  testd <- unique(x$data$testdirection)
  cat("Test direction :", testd, "\n")

  # Convergence info
  if (!is.null(x$fit$convergence) && x$fit$convergence == 0) {
    cat("Optimization : converged\n")
  } else {
    cat("Optimization : NOT converged\n")
  }

  if (!is.null(x$fit$objective)) {
    cat("-2 negative Log Likelihood :", round(x$fit$objective*2, 4), "\n")
  }

  cat("\nUse summary() for parameter estimates.\n")

  invisible(x)
}
