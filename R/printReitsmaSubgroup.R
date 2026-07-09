#' Print Method for ReitsmaSubgroup Objects
#'
#' Displays a concise summary of a fitted Reitsma Subgroup diagnostic test
#' accuracy model, including number of studies, convergence status, and
#' likelihood-based fit statistics.
#'
#' @param x An object of class \code{"ReitsmaSubgroup"}.
#' @param ... Additional arguments (unused).
#'
#' @seealso \code{\link{summary.ReitsmaSubgroup}}
#' @return
#' Invisibly returns the input object.
#' @method print ReitsmaSubgroup
#' @export
print.ReitsmaSubgroup <- function(x, ...) {
  
  cat("\n", "Reitsma Subgroup Model", "\n", sep = "")
  cat(strrep("-", nchar("Reitsma Subgroup Model")), "\n\n", sep = "")
  
  
  n_study   <- nrow(x$data)
  n_sub     <- length(x$subgroups)
  converged <- x$glmmTMB_mu$fit$convergence == 0
  pdHess    <- x$glmmTMB_mu$sdr$pdHess
  ll        <- logLik(x)
  
  cat("Number of studies   :", n_study, "\n")
  cat("Number of subgroups :", n_sub, "\n")
  cat("Optimizer           :", if (converged) "Converged" else "Not converged", "\n")
  cat("Hessian             :", if (pdHess) "Positive definite" else "Not positive definite", "\n")
  cat("Max |grad|          :", max(abs(x$glmmTMB_mu$sdr$gradient.fixed)), "\n")
  cat("-2 log likelihood   :", round(-2 * as.numeric(ll), 3),"( df =", attr(ll, "df"), ")\n")  
  cat("AIC                 :", round(AIC(x), 3), "\n")
  cat("BIC                 :", round(BIC(x), 3), "\n")
  cat("\nUse summary() for parameter estimates.\n")
  
  invisible(x)
}
