#' @keywords internal
#' @noRd
plot_SESPGRID <- function(main){
  plot(1,1, ylim=c(0,1), xlim=c(0,1), xaxt = "n", yaxt="n",
       ann=F, pch=20, col="white",las=1,asp=1)
  axis( side = 1,                          # 1 = bottom axis
        at = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),  # positions of ticks
        labels = c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0))  # custom labels
  axis( side = 2,
        at = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),las=1)
  par(new=TRUE) 
  abline(v=(seq(0,1,0.2)), col="lightgray", lty="dotted")
  abline(h=(seq(0,1,0.2)), col="lightgray", lty="dotted")
  lines(c(0,1),c(0,1),col="lightgray",lty="dotted")
  # Add titles
  title(main=main, xlab="Specificity", ylab="Sensitivity")
}
