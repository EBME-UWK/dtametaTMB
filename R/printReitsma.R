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
  
  cat("\nReitsma Bivariate Random-Effects Model\n")
  cat(strrep("-", 45), "\n")
  
  # Number of studies
  n_study <- nrow(x$data)
  cat("Number of studies:", n_study, "\n")
  
  # Convergence (glmmTMB)
  conv <- tryCatch({
    x$glmmTMB$sdr$pdHess
  }, error = function(e) NA)
  
  conv_msg <- if (isTRUE(conv)) {
    "Converged (positive definite Hessian)"
  } else {
    "Potential convergence issues"
  }
  
  cat("Model fit:", conv_msg, "\n\n")
  
  # Key fixed effects (logit scale)
  est <- x$estimates
  
  cat("Key parameters (logit scale):\n")
  cat("  mu_A (sens):", round(est["mu_A.sens","Estimate"], digits), "\n")
  cat("  mu_B (spec):", round(est["mu_B.spec","Estimate"], digits), "\n")
  
  invisible(x)
}
