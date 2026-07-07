#' Print Method for RutterGatsonis Objects
#'
#' Displays a concise summary of a fitted HSROC model,
#' including number of studies, convergence status, and
#' likelihood-based fit statistics.
#'
#' @param x An object of class \code{"RutterGatsonis"}.
#' @param ... Additional arguments (unused).
#'
#' @seealso \code{\link{summary.RutterGatsonis}}
#' @return
#' Invisibly returns the input object. 
#' @method print RutterGatsonis
#' @export
print.RutterGatsonis <- function(x, ...) {
  
  cat("\n", "Rutter & Gatsonis Model", "\n", sep = "")
  cat(strrep("-", nchar("Rutter & Gatsonis Model")), "\n\n", sep = "")
  
  n_study <- nrow(x$data)
  converged <- x$fit$convergence == 0
  ll <- logLik(x)
  
  cat("Number of studies :", n_study, "\n")
  cat("Model fit         :", if (converged) "Converged" else "Not converged", "\n")
  cat("-2 log likelihood :", round(-2 * as.numeric(ll), 3),"( df =", attr(ll, "df"), ")\n")  
  cat("AIC               :", round(AIC(x), 3), "\n")
  cat("BIC               :", round(BIC(x), 3), "\n")
  cat("\n")
  
  est <- x$sdreport$par.fixed
  
  cat("Lambda            :", round(est["Lambda"], 3), "\n")
  cat("Theta             :", round(est["Theta"], 3), "\n")
  cat("Beta              :", round(est["beta"], 3), "\n")
  
  cat("\nUse summary() for parameter estimates.\n")
  
  invisible(x)
}
