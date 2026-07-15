#' Print Method for RutterGatsonisSubgroupLCA Objects
#'
#' Displays a concise summary of a fitted RutterGatsonisSubgroupLCA model,
#' including number of studies, convergence status, and
#' likelihood-based fit statistics.
#'
#' @param x An object of class \code{"RutterGatsonisSubgroupLCA"}.
#' @param ... Additional arguments (unused).
#'
#' @seealso \code{\link{summary.RutterGatsonisLCA}}
#' @return
#' Invisibly returns the input object. 
#' @method print RutterGatsonisSubgroupLCA
#' @export
print.RutterGatsonisSubgroupLCA <- function(x, ...) {
  
  cat("\n", "Rutter & Gatsonis Subgroup LCA Model", "\n", sep = "")
  cat(strrep("-", nchar("Rutter & Gatsonis Subgroup LCA Model")), "\n\n", sep = "")
  
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