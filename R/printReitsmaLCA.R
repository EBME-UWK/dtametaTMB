#' Print Method for ReitsmaLCA Objects
#'
#' Displays a concise summary of a fitted Reitsma LCA model,
#' including number of studies, convergence status, and
#' likelihood-based fit statistics.
#'
#' @param x An object of class \code{"ReitsmaLCA"}.
#' @param ... Additional arguments (unused).
#'
#' @seealso \code{\link{summary.ReitsmaLCA}}
#' @return
#' Invisibly returns the input object. 
#' @method print ReitsmaLCA
#' @export
print.ReitsmaLCA <- function(x, ...) {
  
  cat("\n", "Reitsma LCA Model", "\n", sep = "")
  cat(strrep("-", nchar("Reitsma LCA Model")), "\n\n", sep = "")
  
  n_study   <- nrow(x$data)
  converged <- x$fit$convergence == 0
  pdHess    <- x$sdreport$pdHess
  ll        <- logLik(x)
  
  cat("Number of studies :", n_study, "\n")
  cat("Optimizer         :", if (converged) "Converged" else "Not converged", "\n")
  cat("Hessian           :", if (pdHess) "Positive definite" else "Not positive definite", "\n")
  cat("Max |grad|        :", max(abs(x$sdreport$gradient.fixed)), "\n")
  cat("-2 log likelihood :", round(-2 * as.numeric(ll), 3),"( df =", attr(ll, "df"), ")\n")  
  cat("AIC               :", round(AIC(x), 3), "\n")
  cat("BIC               :", round(BIC(x), 3), "\n\n")
  
  est <- x$sdreport2[,"Estimate"]
  
  cat("mu_prev           :", round(est["mu_prev"], 3), "\n")
  cat("mu_A.index (sens) :", round(est["mu_A.index"], 3), "\n")
  cat("mu_B.index (spec) :", round(est["mu_B.index"], 3), "\n")
  
  cat("\nUse summary() for parameter estimates.\n")
  
  invisible(x)
}