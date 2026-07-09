#' Print Hoyer AFT Model Object
#'
#' Displays a concise summary of a fitted Hoyer AFT model, including
#' distribution, number of studies, convergence status, and likelihood-based fit statistics.
#'
#' @param x An object of class \code{"HoyerAFT"}.
#' @param ... Further arguments (unused).
#' @return
#' Invisibly returns the input object.
#' @seealso \code{\link{summary.HoyerAFT}}
#' @method print HoyerAFT
#' @export
print.HoyerAFT <- function(x, ...) {
  
  cat("\n", "Hoyer Model", "\n", sep = "")
  cat(strrep("-", nchar("Hoyer Model")), "\n\n", sep = "")
  
  dist_name <- switch(as.character(x$distcode),
                      "1" = "Weibull",
                      "2" = "Lognormal",
                      "3" = "Loglogistic",
                      "Unknown")
  
  n_study   <- length(unique(x$restructured$study))
  converged <- x$fit$convergence == 0
  pdHess    <- x$sdreport$pdHess
  ll        <- logLik(x)

  cat("Number of studies :", n_study, "\n")
  cat("Optimizer         :", if (converged) "Converged" else "Not converged", "\n")
  cat("Hessian           :", if (pdHess) "Positive definite" else "Not positive definite", "\n")
  cat("Max |grad|        :", max(abs(x$sdreport$gradient.fixed)), "\n")
  cat("-2 log likelihood :", round(-2 * as.numeric(ll), 3),"( df =", attr(ll, "df"), ")\n")  
  cat("AIC               :", round(AIC(x), 3), "\n")
  cat("BIC               :", round(BIC(x), 3), "\n")
  cat("\n")
  cat("Distribution      :", dist_name, "\n")
  cat("Test direction    :", unique(x$data$testdirection), "\n")
  
  cat("\nUse summary() for parameter estimates.\n")
  
  invisible(x)
}
