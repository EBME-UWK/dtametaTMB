#' Print a diagnostic test accuracy model summary
#'
#' Prints summary objects returned by
#' \code{summary()} methods for diagnostic
#' test accuracy models.
#'
#' @param x A summary object inheriting from
#'   \code{"summary.DTAmodel"}.
#' @param digits Number of significant digits used for printing.
#' @param ... Additional arguments passed to printing methods.
#'
#' @return Invisibly returns \code{x}.
#' @export

print.summary.DTAmodel <- function(x,
                                   digits = max(8L, getOption("digits") - 8L),
                                   ...) {
  
  mod <- c(
    summary.HoyerAFT = "Summary: Hoyer Model",
    summary.Reitsma = "Summary: Reitsma Model",
    summary.ReitsmaLCA = "Summary: Reitsma LCA Model",
    summary.ReitsmaSubgroup = "Summary: Reitsma Subgroup Model",
    summary.ReitsmaSubgroupLCA = "Summary: Reitsma Subgroup LCA Model",
    summary.RutterGatsonis = "Summary: Rutter & Gatsonis Model",
    summary.RutterGatsonisReg = "Summary: Rutter & Gatsonis Regression Model",
    summary.RutterGatsonisLCA = "Summary: Rutter & Gatsonis LCA Model",
    summary.RutterGatsonisSubgroup = "Summary: Rutter & Gatsonis Subgroup Model",
    summary.RutterGatsonisSubgroupLCA = "Summary: Rutter & Gatsonis Subgroup LCA Model"
  )
  
  titles <- c(
    estimates                 = "Parameter estimates",
    sdreport2                 = "Parameter estimates",
    sensspec                  = "Sensitivity / specificity",
    RutterGatsonis_recovered  = "Recovered HSROC parameters (Rutter-Gatsonis)",
    Reitsma_recovered         = "Recovered bivariate parameters (Reitsma)",
    prevref                   = "Prevalence and reference standard",
    subgroups                 = "Subgroups"
  )
  
  cl     <- class(x)[1]
  header <- if (cl %in% names(mod)) mod[[cl]] else paste0("Summary: ", sub("^summary\\.", "", cl), " Model")
  
  cat("\n", header, "\n", sep = "")
  cat(strrep("-", nchar(header)), "\n", sep = "")
  
  for (nm in names(x)) {
    
    val <- x[[nm]]
    if (is.null(val)) next
    
    title <- if (nm %in% names(titles)) titles[[nm]] else getdtaSectionTitle(nm)
    
    cat("\n", title, "\n", sep = "")
    cat(strrep("-", nchar(title)), "\n", sep = "")
    
    if (is.data.frame(val) || is.matrix(val)) {
      print(getdtaRoundNumeric(val, digits), digits = digits)
    } else {
      print(val)
    }
  }
  
  cat("\n")
  invisible(x)
}

#' @keywords internal
#' @noRd
getdtaSectionTitle <- function(nm) {
  nm <- gsub("[._]+", " ", nm)
  nm <- gsub("([a-z])([A-Z])", "\\1 \\2", nm)
  nm <- trimws(nm)
  paste0(toupper(substring(nm, 1, 1)), substring(nm, 2))
}

#' @keywords internal
#' @noRd
getdtaRoundNumeric <- function(val, digits) {
  if (is.matrix(val)) {
    if (is.numeric(val)) return(signif(val, digits))
    return(val)
  }
  for (j in seq_along(val)) {
    if (is.numeric(val[[j]])) val[[j]] <- signif(val[[j]], digits)
  }
  val
}