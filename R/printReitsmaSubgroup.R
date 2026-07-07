#' Print Method for ReitsmaSubgroup Objects
#'
#' Displays a concise summary of a fitted Reitsma Subgroup diagnostic test
#' accuracy model, including number of studies, convergence status.
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
  
  
  n_study <- nrow(x$data)
  n_sub   <- length(x$subgroups)
  
  converged <- tryCatch({
    isTRUE(x$glmmTMB_mu$sdr$pdHess)
  }, error = function(e) FALSE)
  ll <- logLik(x)
  
  cat("Number of studies   :", n_study, "\n")
  cat("Number of subgroups :", n_sub, "\n")
  cat("Model fit           :", if (converged) "Converged" else "Not converged", "\n")
  cat("-2 log likelihood   :", round(-2 * as.numeric(ll), 3),"( df =", attr(ll, "df"), ")\n")  
  cat("AIC                 :", round(AIC(x), 3), "\n")
  cat("BIC                 :", round(BIC(x), 3), "\n")
  cat("\n")
  
  cat("\nUse summary() for parameter estimates.\n")
  
  invisible(x)
}
