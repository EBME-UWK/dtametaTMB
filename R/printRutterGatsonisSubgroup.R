#' Print Method for RutterGatsonisSubgroup Objects
#'
#' Displays a concise summary of a fitted HSROC model,
#' including number of studies, and convergence status.
#'
#' @param x An object of class \code{"RutterGatsonisSubgroup"}.
#' @param ... Additional arguments (unused).
#'
#' @seealso \code{\link{summary.RutterGatsonisSubgroup}}
#'
#' @method print RutterGatsonisSubgroup
#' @export
print.RutterGatsonisSubgroup <- function(x, ...) {
  
  cat("\n", "Rutter & Gatsonis Subgroup Model", "\n", sep = "")
  cat(strrep("-", nchar("Rutter & Gatsonis Subgroup Model")), "\n\n", sep = "")
  
  n_study <- nrow(x$data)
  n_sub   <- length(x$subgroups)
  converged <- x$fit$convergence == 0
  loglik <- if (!is.null(x$fit$objective)) 2 * x$fit$objective else NULL
  
  
  cat("Number of studies   :", n_study, "\n")
  cat("Number of subgroups :", n_sub, "\n")
  cat("Model fit           :", if (converged) "Converged" else "Not converged", "\n")
  
  if (!is.null(loglik)) {
  cat("-2 log likelihood   :", round(loglik, 3), "\n")
  }
  cat("AIC                 :", round(AIC(x), 3), "\n")
  cat("BIC                 :", round(BIC(x), 3), "\n")
  cat("\n")
  
  cat("\nUse summary() for parameter estimates.\n")
  
  invisible(x)
}
