#' Print Method for ReitsmaSubgroupLCA Objects
#'
#' Displays a concise summary of a fitted ReitsmaSubgroupLCA model,
#' including number of studies, convergence status, and
#' likelihood-based fit statistics.
#'
#' @param x An object of class \code{"ReitsmaSubgroupLCA"}.
#' @param ... Additional arguments (unused).
#'
#' @seealso \code{\link{summary.ReitsmaLCA}}
#' @return
#' Invisibly returns the input object. 
#' @method print ReitsmaSubgroupLCA
#' @export
print.ReitsmaSubgroupLCA <- function(x, ...) {
  
  cat("\n", "Reitsma Subgroup LCA Model", "\n", sep = "")
  cat(strrep("-", nchar("Reitsma Subgroup LCA Model")), "\n\n", sep = "")
  
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
  cat("BIC               :", round(BIC(x), 3), "\n")
  
  cat("\nUse summary() for parameter estimates.\n")
  
  invisible(x)
}