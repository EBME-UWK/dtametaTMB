#' Print Method for Reitsma Objects
#'
#' Displays a concise summary of a fitted Reitsma diagnostic test
#' accuracy model, including number of studies, convergence status
#' and likelihood-based fit statistics.
#'
#' @param x An object of class \code{"Reitsma"}.
#' @param ... Additional arguments (unused).
#'
#' @seealso \code{\link{summary.Reitsma}}
#' @return
#' Invisibly returns the input object.
#' @method print Reitsma
#' @export
print.Reitsma <- function(x, ...) {
  
  cat("\n", "Reitsma Model", "\n", sep = "")
  cat(strrep("-", nchar("Reitsma Model")), "\n\n", sep = "")
  
  
  n_study   <- nrow(x$data)
  converged <- x$glmmTMB$fit$convergence == 0
  pdHess    <- x$glmmTMB$sdr$pdHess
  ll        <- logLik(x)
  
  cat("Number of studies :", n_study, "\n")
  cat("Optimizer         :", if (converged) "Converged" else "Not converged", "\n")
  cat("Hessian           :", if (pdHess) "Positive definite" else "Not positive definite", "\n")
  cat("Max |grad|        :", max(abs(x$glmmTMB$sdr$gradient.fixed)), "\n")
  cat("-2 log likelihood :", round(-2 * as.numeric(ll), 3),"( df =", attr(ll, "df"), ")\n")  
  cat("AIC               :", round(AIC(x), 3), "\n")
  cat("BIC               :", round(BIC(x), 3), "\n")
  
  cat("\n")
  
  est <- x$estimates
  
  cat("mu_A (logitsens)  :", round(est["mu_A.sens","Estimate"], 3), "\n")
  cat("mu_B (logitspec)  :", round(est["mu_B.spec","Estimate"], 3), "\n")
  
  cat("\nUse summary() for parameter estimates.\n")

  invisible(x)
}
