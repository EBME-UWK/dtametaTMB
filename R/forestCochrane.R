#' Coupled Forest plot for diagnostic test accuracy meta-analysis
#'
#' @param x Object of class \code{"Cochrane"} such as \code{"RutterGatsonis"}, \code{"Reitsma"} or \code{"HoyerAFT"}
#' @param ... Additional graphical arguments (not currently in use)
#'
#' @method forest Cochrane
#' @importFrom forestploter forest
#' @importFrom grid unit 
#' @importFrom stats qbeta
#' @export
forest.Cochrane <- function(x, ...) {
  XP <- x$data
  XP$FPR  <- 1 - XP$spec
  # How do I get Clopper-Pearson confidence limits for sensitivity and specificities?
  XP$Sens_LCI  <- with(XP,stats::qbeta(p=0.05/2,  shape1=TP,  shape2=FN+1))
  XP$Sens_UCI  <- with(XP,stats::qbeta(p=1-0.05/2,shape1=TP+1,shape2=FN))
  XP$FPR_LCI   <- with(XP,stats::qbeta(p=0.05/2,  shape1=FP,  shape2=TN+1))
  XP$FPR_UCI   <- with(XP,stats::qbeta(p=1-0.05/2,shape1=FP+1,shape2=TN))
  # How do I create the forest plot?
  XP$`Sensitivity (95%-CI)` <- with(XP,paste0(sprintf("%.2f", sens)," [",
                                              sprintf("%.2f", Sens_LCI),", ",
                                              sprintf("%.2f", Sens_UCI),"]"))
  XP$`Specificity (95%-CI)` <- with(XP,paste0(sprintf("%.2f", spec)," [",
                                              sprintf("%.2f", 1-FPR_UCI),", ",
                                              sprintf("%.2f", 1-FPR_LCI),"]"))
  
  if (inherits(x, "Reitsma") || inherits(x, "RutterGatsonis")) {
    dt <- XP[,c("study","TP","FP","FN","TN","Sensitivity (95%-CI)","Specificity (95%-CI)")]
    dt$" "    <- " "
    dt$fsens  <- paste(rep(" ",18),collapse=" ")
    dt$a      <- " "
    dt$fspec  <- paste(rep(" ",18),collapse=" ")  
    cc <- colnames(dt) 
    colnames(dt) <- c("Study",cc[2:8],"Sensitivity (95%-CI)"," ","Specificity (95%-CI)")
    
    p <- forestploter::forest(dt,
                              est = list(XP$sens,
                                         XP$spec),
                              lower = list(XP$Sens_LCI,
                                           1-XP$FPR_UCI), 
                              upper = list(XP$Sens_UCI,
                                           1-XP$FPR_LCI),
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
  }
  if(inherits(x,"HoyerAFT")){
    dt <- XP[, c("study", "threshold", "TP","FP","FN","TN", "Sensitivity (95%-CI)", "Specificity (95%-CI)")]
    
    dt <- dt[order(dt$study, dt$threshold), ]
    dt$study <- as.character(dt$study)
    dt$study[duplicated(dt$study)] <- paste0("   ", dt$Study[duplicated(dt$study)])
    
    dt$" "    <- " "
    dt$fsens  <- paste(rep(" ",18),collapse=" ")
    dt$a      <- " "
    dt$fspec  <- paste(rep(" ",18),collapse=" ")  
    cc <- colnames(dt) 
    colnames(dt) <- c("Study","Threshold",cc[3:9],"Sensitivity (95%-CI)"," ","Specificity (95%-CI)")
    
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
                                 col = 2:8,
                                 which="text",
                                 hjust = grid::unit(1,"npc"),
                                 x = grid::unit(1,"npc"))
    p <- forestploter::edit_plot(p,
                                 col = 2:12,
                                 part="header",
                                 hjust = grid::unit(1,"npc"),
                                 x = grid::unit(1,"npc"))
    plot(p)
    
  }
  

}
  
  