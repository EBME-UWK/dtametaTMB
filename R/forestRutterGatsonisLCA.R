#' Coupled Forest plot for latent class diagnostic test accuracy meta-analysis
#'
#' Study-specific sensitivities and specificities are empirical Bayes estimates
#' derived from the fitted latent class model. Confidence intervals are based on 
#' the conditional posterior variance of the study-specific random effects 
#' obtained from the TMB Laplace approximation.
#' 
#' @param x Object of class \code{"RutterGatsonisLCA"}
#' @param conflevel Confidence level for confidence intervals. Default is 0.95.
#' @param ... Additional graphical arguments (not currently in use)
#'
#' @method forest RutterGatsonisLCA
#' @importFrom forestploter forest edit_plot
#' @importFrom grid unit
#' @importFrom stats qnorm plogis qlogis
#' @return
#' Invisibly returns a \code{forestploter} object. Users may further modify the plot
#' using \code{forestploter} functions before printing or exporting.
#' @export
forest.RutterGatsonisLCA <- function(x,conflevel=0.95, ...) {
  if (!is.numeric(conflevel) || length(conflevel) != 1L ||
      conflevel <= 0 || conflevel >= 1) {
    stop("conflevel must be a single number in (0, 1).")
  }
  XP <- x$data
  alpha   <- 1 - conflevel
  qq      <- stats::qnorm(1-alpha/2)
  senslabel <-  paste0("Sensitivity (",round(100 * conflevel), "%-CI)")
  speclabel <-  paste0("Specificity (",round(100 * conflevel), "%-CI)")
  # How do I get confidence limits for sensitivity and specificities?
  lsens        <- stats::qlogis(XP$sens_eb)
  lspec        <- stats::qlogis(XP$spec_eb)
  XP$Sens_LCI  <- with(XP,stats::plogis(lsens-qq*sqrt(lsens_eb_var)))
  XP$Sens_UCI  <- with(XP,stats::plogis(lsens+qq*sqrt(lsens_eb_var)))
  XP$Spec_LCI  <- with(XP,stats::plogis(lspec-qq*sqrt(lspec_eb_var)))
  XP$Spec_UCI  <- with(XP,stats::plogis(lspec+qq*sqrt(lspec_eb_var)))
  # How do I create the forest plot?
  XP$senslabel <- with(XP,paste0(sprintf("%.2f", sens_eb)," [",
                                 sprintf("%.2f", Sens_LCI),", ",
                                 sprintf("%.2f", Sens_UCI),"]"))
  XP$speclabel <- with(XP,paste0(sprintf("%.2f", spec_eb)," [",
                                 sprintf("%.2f", Spec_LCI),", ",
                                 sprintf("%.2f", Spec_UCI),"]"))
  
  XP <- XP[order(XP$study), ]
  dt <- XP[,c("study","y11","y10","y01","y00","senslabel","speclabel")]
  dt$" "    <- " "
  dt$fsens  <- paste(rep(" ",18),collapse=" ")
  dt$a      <- " "
  dt$fspec  <- paste(rep(" ",18),collapse=" ")  
  cc <- colnames(dt) 
  colnames(dt) <- c("Study","+/+","+/-","-/+","-/-",senslabel,speclabel," ",senslabel," ",speclabel)
    
  p <- forestploter::forest(dt,
                            est = list(XP$sens_eb,
                                       XP$spec_eb),
                            lower = list(XP$Sens_LCI,
                                         XP$Spec_LCI), 
                            upper = list(XP$Sens_UCI,
                                         XP$Spec_UCI),
                            sizes = 0.75,
                            ci_column = c(9,11),
                            nudge_y=0.000001,
                            xlim=c(0,1),
                            ref_line = 3)
  p <- forestploter::edit_plot(p,
                               col = 2:7,
                               which="text",
                               hjust = grid::unit(1,"npc"),
                               x = grid::unit(1,"npc"))
  p <- forestploter::edit_plot(p,
                               col = 2:11,
                               part="header",
                               hjust = grid::unit(1,"npc"),
                               x = grid::unit(1,"npc"))
  plot(p)
  invisible(p)
}

