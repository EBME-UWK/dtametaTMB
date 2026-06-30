#' Print Method for RutterGatsonis Objects
#'
#' Displays a concise summary of a fitted HSROC model,
#' including number of studies, convergence status, and key estimates.
#'
#' @param x An object of class \code{"RutterGatsonis"}.
#' @param ... Additional arguments (unused).
#'
#' @seealso \code{\link{summary.RutterGatsonis}}
#'
#' @method print RutterGatsonis
#' @export
print.RutterGatsonis <- function(x, ...) {
  
  cat("\n", "Rutter & Gatsonis Model", "\n", sep = "")
  cat(strrep("-", nchar("Rutter & Gatsonis Model")), "\n\n", sep = "")
  
  n_study <- nrow(x$data)
  converged <- x$fit$convergence == 0
  loglik <- if (!is.null(x$fit$objective)) 2 * x$fit$objective else NULL
  
  
  cat("Number of studies :", n_study, "\n")
  cat("Model fit         :", if (converged) "Converged" else "Not converged", "\n")
  
  if (!is.null(loglik)) {
  cat("-2 log likelihood :", round(loglik, 3), "\n")
  }
  cat("\n")
  
  est <- x$sdreport$par.fixed
  
  cat("Lambda            :", round(est["Lambda"], 3), "\n")
  cat("Theta             :", round(est["Theta"], 3), "\n")
  cat("Beta              :", round(est["beta"], 3), "\n")
  
  cat("\nUse summary() for parameter estimates.\n")
  
  invisible(x)
}
