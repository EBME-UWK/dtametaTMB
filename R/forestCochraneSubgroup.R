#' Coupled Forest plot for diagnostic test accuracy meta-analysis
#' 
#' Provides coupled forest plots of sensitivities and specificities
#' with Clopper-Pearson confidence limits.
#'
#' @param x Object of class \code{"CochraneSubgroup"} such as \code{"RutterGatsonisSubgroup"}, \code{"ReitsmaSubgroup"}
#' @param conflevel Confidence level for confidence intervals. Default is 0.95.
#' @param subgroup_label Column name for the subgroup. Defaults to \code{"Subgroup"}.
#' @param ... Additional graphical arguments (not currently in use)
#'
#' @method forest CochraneSubgroup
#' @importFrom forestploter forest
#' @importFrom grid unit 
#' @importFrom stats qbeta
#' @return
#' No return value. Called for its side effect of producing a plot.
#' @export
forest.CochraneSubgroup <- function(x, conflevel=0.95, subgroup_label="Subgroup", ...) {
  if (!is.numeric(conflevel) || length(conflevel) != 1L ||
      conflevel <= 0 || conflevel >= 1) {
    stop("conflevel must be a single number in (0, 1).")
  }
  XP <- x$data
  XP$FPR  <- 1 - XP$spec
  alpha   <- 1 - conflevel
  senslabel <-  paste0("Sensitivity (",round(100 * conflevel), "%-CI)")
  speclabel <-  paste0("Specificity (",round(100 * conflevel), "%-CI)")
  # How do I get Clopper-Pearson confidence limits for sensitivity and specificities?
  XP$Sens_LCI  <- with(XP,stats::qbeta(p=alpha/2,  shape1=TP,  shape2=FN+1))
  XP$Sens_UCI  <- with(XP,stats::qbeta(p=1-alpha/2,shape1=TP+1,shape2=FN))
  XP$FPR_LCI   <- with(XP,stats::qbeta(p=alpha/2,  shape1=FP,  shape2=TN+1))
  XP$FPR_UCI   <- with(XP,stats::qbeta(p=1-alpha/2,shape1=FP+1,shape2=TN))
  # How do I create the forest plot?
  XP$senslabel <- with(XP,paste0(sprintf("%.2f", sens)," [",
                                         sprintf("%.2f", Sens_LCI),", ",
                                         sprintf("%.2f", Sens_UCI),"]"))
  XP$speclabel <- with(XP,paste0(sprintf("%.2f", spec)," [",
                                         sprintf("%.2f", 1-FPR_UCI),", ",
                                         sprintf("%.2f", 1-FPR_LCI),"]"))
  
  if (inherits(x, "ReitsmaSubgroup") || inherits(x, "RutterGatsonisSubgroup")) {
    XP <- XP[order(XP$study,XP$subgroup), ]
    dt <- XP[,c("study","subgroup","TP","FP","FN","TN","senslabel","speclabel")]
    dt$" "    <- " "
    dt$fsens  <- paste(rep(" ",18),collapse=" ")
    dt$a      <- " "
    dt$fspec  <- paste(rep(" ",18),collapse=" ")  
    cc <- colnames(dt) 
    colnames(dt) <- c("Study",subgroup_label,cc[3:6],
                      senslabel,speclabel," ",senslabel," ",speclabel)
    
    p <- forestploter::forest(dt,
                              est = list(XP$sens,
                                         XP$spec),
                              lower = list(XP$Sens_LCI,
                                           1-XP$FPR_UCI), 
                              upper = list(XP$Sens_UCI,
                                           1-XP$FPR_LCI),
                              sizes = 0.75,
                              ci_column = c(10,12),
                              nudge_y=0.000001,
                              xlim=c(0,1),
                              ref_line = 3)
    p <- forestploter::edit_plot(p,
                                 col = 3:8,
                                 which="text",
                                 hjust = grid::unit(1,"npc"),
                                 x = grid::unit(1,"npc"))
    p <- forestploter::edit_plot(p,
                                 col = 3:12,
                                 part="header",
                                 hjust = grid::unit(1,"npc"),
                                 x = grid::unit(1,"npc"))
    plot(p)
  }
  # if(inherits(x,"HoyerAFT")){
  #   XP <- XP[order(XP$study,XP$subgroup,XP$threshold), ]
  #   dt <- XP[, c("study",subgroup_label,"threshold", "TP","FP","FN","TN", "senslabel", "speclabel")]
  #   
  #   dt$study <- as.character(dt$study)
  #   dt$study[duplicated(dt$study)] <- " "
  #   dt$subgroup[duplicated(dt$study)] <- " "
  #   dt$study[duplicated(dt$study)] <- " "
  #   
  #   dt$" "    <- " "
  #   dt$fsens  <- paste(rep(" ",18),collapse=" ")
  #   dt$a      <- " "
  #   dt$fspec  <- paste(rep(" ",18),collapse=" ")  
  #   cc <- colnames(dt) 
  #   colnames(dt) <- c("Study","Subgroup","Threshold",cc[4:7],
  #                     senslabel,speclabel," ", senslabel," ",speclabel)
  #   
  #   p <- forestploter::forest(dt,
  #                             est = list(XP$sens,
  #                                        XP$spec),
  #                             lower = list(XP$Sens_LCI,
  #                                          1-XP$FPR_UCI), 
  #                             upper = list(XP$Sens_UCI,
  #                                          1-XP$FPR_LCI),
  #                             sizes = 0.75,
  #                             ci_column = c(11,13),
  #                             nudge_y=0.000001,
  #                             xlim=c(0,1),
  #                             ref_line = 3)
  #   p <- forestploter::edit_plot(p,
  #                                col = 2:9,
  #                                which="text",
  #                                hjust = grid::unit(1,"npc"),
  #                                x = grid::unit(1,"npc"))
  #   p <- forestploter::edit_plot(p,
  #                                col = 2:13,
  #                                part="header",
  #                                hjust = grid::unit(1,"npc"),
  #                                x = grid::unit(1,"npc"))
  #   plot(p)
  #   
  # }
  invisible(NULL)
}

