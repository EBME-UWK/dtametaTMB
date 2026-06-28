#' Print Method for Reitsma Objects
#'
#' Displays a concise summary of a fitted Reitsma diagnostic test
#' accuracy model, including number of studies, convergence status,
#' and key parameter estimates.
#'
#' @param x An object of class \code{"Reitsma"}.
#' @param digits Number of digits to print (default: 3).
#' @param ... Additional arguments (unused).
#'
#' @return
#' Invisibly returns the input object.
#'
#' @seealso \code{\link{summary.Reitsma}}
#'
#' @method print Reitsma
#' @export
print.Reitsma <- function(x, digits = 3, ...) {
  
  cat("\n", "Reitsma Model", "\n", sep = "")
  cat(strrep("-", nchar("Reitsma Model")), "\n\n", sep = "")
  
  
  n_study <- nrow(x$data)
  
  converged <- tryCatch({
    isTRUE(x$glmmTMB$sdr$pdHess)
  }, error = function(e) FALSE)
  
  # logLik extraction (if available)
  loglik <- tryCatch({
    -2 * as.numeric(stats::logLik(x$glmmTMB))
  }, error = function(e) NULL)
  
  
  cat("Number of studies :", n_study, "\n")
  cat("Model fit         :", if (converged) "Converged" else "Not converged", "\n")
  if (!is.null(loglik)) {
    cat("-2 log likelihood :", round(loglik, 3), "\n")
  }
  cat("\n")
  
  est <- x$estimates
  
  cat("mu_A (sens)       :", round(est["mu_A.sens","Estimate"], 3), "\n")
  cat("mu_B (spec)       :", round(est["mu_B.spec","Estimate"], 3), "\n")
  
  cat("\nUse summary() for parameter estimates.\n")

  invisible(x)
}
