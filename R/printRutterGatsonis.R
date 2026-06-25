#' Print Method for RutterGatsonis Objects
#'
#' Displays a concise summary of a fitted HSROC model,
#' including number of studies, convergence status, and key estimates.
#'
#' @param x An object of class \code{"RutterGatsonis"}.
#' @param digits Number of digits to print (default: 3).
#' @param ... Additional arguments (unused).
#'
#' @return
#' Invisibly returns the input object.
#'
#' @seealso \code{\link{summary.RutterGatsonis}}
#'
#' @method print RutterGatsonis
#' @export
print.RutterGatsonis <- function(x, digits = 3, ...) {
  
  cat("\nRutter & Gatsonis HSROC Model\n")
  cat(strrep("-", 40), "\n")
  
  # Number of studies
  n_study <- nrow(x$data)
  cat("Number of studies:", n_study, "\n")
  
  # Convergence
  conv <- x$fit$convergence
  conv_msg <- if (conv == 0) "Converged" else paste("Not converged (code:", conv, ")")
  cat("Model fit:", conv_msg, "\n\n")
  
  # Key parameters
  est <- x$sdreport$par.fixed
  
  cat("Parameters (HSROC scale):\n")
  cat("  Lambda:", round(est["Lambda"], digits), "\n")
  cat("  Theta :", round(est["Theta"], digits), "\n")
  cat("  Beta  :", round(est["beta"], digits), "\n")
  
  invisible(x)
}
