#' Print Method for RutterGatsonisReg Objects
#'
#' Displays a concise summary of a fitted HSROC model,
#' including number of studies, convergence status,
#' and likelihood-based fit statistics.
#'
#' @param x An object of class \code{"RutterGatsonisReg"}.
#' @param ... Additional arguments (unused).
#'
#' @seealso \code{\link{summary.RutterGatsonisReg}}
#' @return Invisibly returns the input object.
#' @method print RutterGatsonisReg
#' @export
print.RutterGatsonisReg <- function(x, ...) {
  
  cat("\n", "Rutter & Gatsonis Regression Model", "\n", sep = "")
  cat(strrep("-", nchar("Rutter & Gatsonis Regression Model")), "\n\n", sep = "")
  
  n_study <- nrow(x$data)
  converged <- x$fit$convergence == 0
  ll <- logLik(x)
  
  cat("Number of studies :", n_study, "\n")
  cat("Model fit         :", if (converged) "Converged" else "Not converged", "\n")
  cat("-2 log likelihood :", round(-2 * as.numeric(ll), 3),"( df =", attr(ll, "df"), ")\n")  
  cat("AIC               :", round(AIC(x), 3), "\n")
  cat("BIC               :", round(BIC(x), 3), "\n")
  cat("\n")
  
  cat("\nUse summary() for parameter estimates.\n")
  
  invisible(x)
}
