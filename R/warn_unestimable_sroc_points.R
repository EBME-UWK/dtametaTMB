#' @keywords internal
#' @noRd
warn_unestimable_sroc_points <- function(data) {
  no_diseased <- data$n1 == 0
  no_nondiseased <- data$n0 == 0
  
  affected <- no_diseased | no_nondiseased
  
  if(!any(affected)) {
    return(invisible(NULL))
  }
  
  reasons <- character(sum(affected))
  reasons[no_diseased[affected] & !no_nondiseased[affected]] <- "no diseased participants"
  reasons[!no_diseased[affected] & no_nondiseased[affected]] <- "no non-diseased participants"
  reasons[no_diseased[affected] & no_nondiseased[affected]]  <- "neither diseased nor non-diseased participants"
  details <- paste0( "\"", data$study[affected], "\""," (", reasons, ")")
  
  warning(
    "The following studies are not displayed in the SROC plot ",
    "because sensitivity or specificity is not estimable: ",
    paste(details, collapse = "; "),
    "."
  )
  invisible(affected)
}