#' @keywords internal
#' @noRd

plotStudyLevelEstimates <- function(xsens,
                                    xspec,
                                    xdata,
                                    symbol,
                                    pctse,
                                    pctsp,
                                    scale,
                                    col,
                                    xstud,
                                    conflevel,
                                    studyCI,
                                    studylabels,
                                    LCA=FALSE){
  
   pointsXY(x = 1-xspec, 
            y = xsens, 
            symbol = symbol, 
            scale = scale*0.5,
            cex.x = pctsp,
            cex.y = pctse,
            col=col)
   symb <- switch(symbol, rectangle = 0, plus = 3, cross = 4, star = 8, ellipse = 1, diamond = 5, triangle = 2)
   if(isTRUE(studyCI)) {
     if(!isTRUE(LCA)){
       forestci <- getForestSensSpec(xdata=xdata,
                                     conflevel=conflevel)$XP[,c("sens",
                                                                "spec",
                                                                "Sens_LCI",
                                                                "Sens_UCI",
                                                                "Spec_LCI",
                                                                "Spec_UCI")]
     }
     if(isTRUE(LCA)){
       forestci <- getForestSensSpecLCA(xdata=xdata,
                                        conflevel=conflevel)$XP[,c("sens",
                                                                   "spec",
                                                                   "Sens_LCI",
                                                                   "Sens_UCI",
                                                                   "Spec_LCI",
                                                                   "Spec_UCI")]
     }

     drawStudyPointCI(x = 1 - forestci$spec,
                      y = forestci$sens,
                      sens_lower = forestci$Sens_LCI,
                      sens_upper = forestci$Sens_UCI,
                      spec_lower = forestci$Spec_LCI,
                      spec_upper = forestci$Spec_UCI,
                      symbol = symbol,
                      col = col,
                      scale = scale * 0.5,
                      cex.x = pctsp,
                      cex.y = pctse)
   }
   if(isTRUE(studylabels)){
      graphics::text(x = 1-xspec,
                     y = xsens,
                     cex = 1,
                     pos = 4,
                     col = col,
                     labels = xstud)
   }
   return(symb)
}