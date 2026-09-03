#' Forest plot generic
#'
#' Produces coupled forest plots.
#'
#' @param x Object
#' @param ... Additional arguments
#'
#' @return
#' Invisibly returns a \code{forestploter} object. Users may further modify the plot
#' using \code{forestploter} functions before printing or exporting.
#' @export
forest <- function(x, ...) {
  UseMethod("forest")
}

#' @keywords internal
#' @importFrom stats qbeta
#' @noRd
getForestSensSpec <- function(x,conflevel){ 
  XP <- x$data
  XP$FPR  <- 1 - XP$spec
  alpha   <- 1 - conflevel
  senslab <-  paste0("Sensitivity (",round(100 * conflevel), "%-CI)")
  speclab <-  paste0("Specificity (",round(100 * conflevel), "%-CI)")
  # How do I get Clopper-Pearson confidence limits for sensitivity and specificities?
  XP$Sens_LCI  <- with(XP,stats::qbeta(p=alpha/2,  shape1=TP,  shape2=FN+1))
  XP$Sens_UCI  <- with(XP,stats::qbeta(p=1-alpha/2,shape1=TP+1,shape2=FN))
  XP$Spec_UCI  <- with(XP,1-stats::qbeta(p=alpha/2,  shape1=FP,  shape2=TN+1))
  XP$Spec_LCI  <- with(XP,1-stats::qbeta(p=1-alpha/2,shape1=FP+1,shape2=TN))
  # How do I create the forest plot?
  XP$senslabel <- with(XP,paste0(sprintf("%.2f", sens)," [",
                                 sprintf("%.2f", Sens_LCI),", ",
                                 sprintf("%.2f", Sens_UCI),"]"))
  XP$speclabel <- with(XP,paste0(sprintf("%.2f", spec)," [",
                                 sprintf("%.2f", Spec_LCI),", ",
                                 sprintf("%.2f", Spec_UCI),"]"))
  return(list(XP=XP,
              senslab=senslab,
              speclab=speclab))
}

#' @keywords internal
#' @importFrom stats qnorm qlogis plogis
#' @noRd
getForestSensSpecLCA <- function(x,conflevel){ 
  XP <- x$data
  alpha   <- 1 - conflevel
  qq      <- stats::qnorm(1-alpha/2)
  senslab <-  paste0("Sensitivity (",round(100 * conflevel), "%-CI)")
  speclab <-  paste0("Specificity (",round(100 * conflevel), "%-CI)")
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
  XP$sens <- XP$sens_eb
  XP$spec <- XP$spec_eb
  return(list(XP=XP,
              senslab=senslab,
              speclab=speclab))
}


#' @keywords internal
#' @importFrom forestploter forest edit_plot
#' @importFrom grid unit
#' @noRd

getForestPlot <- function(dt,XP){
  p <- forestploter::forest(dt,
                            est = list(XP$sens,
                                       XP$spec),
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

#' @keywords internal
#' @importFrom forestploter forest edit_plot
#' @importFrom grid unit
#' @noRd

getForestPlotSub <- function(dt,XP,s=3){
  p <- forestploter::forest(dt,
                            est = list(XP$sens,
                                       XP$spec),
                            lower = list(XP$Sens_LCI,
                                         XP$Spec_LCI), 
                            upper = list(XP$Sens_UCI,
                                         XP$Spec_UCI),
                            sizes = 0.75,
                            ci_column = c(10,12),
                            nudge_y=0.000001,
                            xlim=c(0,1),
                            ref_line = 3)
  p <- forestploter::edit_plot(p,
                               col = s:8,
                               which="text",
                               hjust = grid::unit(1,"npc"),
                               x = grid::unit(1,"npc"))
  p <- forestploter::edit_plot(p,
                               col = s:12,
                               part="header",
                               hjust = grid::unit(1,"npc"),
                               x = grid::unit(1,"npc"))
  plot(p)
  invisible(p)
}