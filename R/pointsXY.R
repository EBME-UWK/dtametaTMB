#' @keywords internal
#' @importFrom graphics polygon segments
#' @noRd
pointsXY <- function(x, y,
                     symbol = "rectangle",
                     scale = 1,
                     cex.x = 1,
                     cex.y = 1,
                     col = par("fg"),
                     lwd = 1) {
  
  stopifnot(length(x) == length(y) & 
            length(x) == length(cex.x) &
            length(x) == length(cex.y))
  
  for (i in seq_along(x)) {
    
    rx <- cex.x[i] * scale
    ry <- cex.y[i] * scale
    
    ## PLUS ---------------------------------------------------------
    if (symbol == "plus") {
      
      graphics::segments(
        x[i] - rx, y[i],
        x[i] + rx, y[i],
        col = col, lwd = lwd
      )
      
      graphics::segments(
        x[i], y[i] - ry,
        x[i], y[i] + ry,
        col = col, lwd = lwd
      )
      
      ## CROSS --------------------------------------------------------
    } else if (symbol == "cross") {
      
      graphics::segments(
        x[i] - rx, y[i] - ry,
        x[i] + rx, y[i] + ry,
        col = col, lwd = lwd
      )
      
      graphics::segments(
        x[i] - rx, y[i] + ry,
        x[i] + rx, y[i] - ry,
        col = col, lwd = lwd
      )
      
      ## STAR = PLUS + CROSS -----------------------------------------
    } else if (symbol == "star") {
      
      # horizontal
      graphics::segments(
        x[i] - rx, y[i],
        x[i] + rx, y[i],
        col = col, lwd = lwd
      )
      
      # vertical
      graphics::segments(
        x[i], y[i] - ry,
        x[i], y[i] + ry,
        col = col, lwd = lwd
      )
      
      # diagonal /
      graphics::segments(
        x[i] - rx, y[i] - ry,
        x[i] + rx, y[i] + ry,
        col = col, lwd = lwd
      )
      
      # diagonal \
      graphics::segments(
        x[i] - rx, y[i] + ry,
        x[i] + rx, y[i] - ry,
        col = col, lwd = lwd
      )
      
      ## ELLIPSE -----------------------------------------------------
    } else if (symbol == "ellipse") {
      
      theta <- seq(0, 2*pi, length.out = 100)
      
      graphics::polygon(
        x[i] + rx * cos(theta),
        y[i] + ry * sin(theta),
        col = NA,
        border = col,
        lwd = lwd
      )
      
      ## SQUARE -------------------------------------------------------
    } else if (symbol == "rectangle") {
      
      graphics::polygon(
        c(x[i] - rx, x[i] + rx,
          x[i] + rx, x[i] - rx),
        c(y[i] - ry, y[i] - ry,
          y[i] + ry, y[i] + ry),
        col = NA,
        border = col,
        lwd = lwd
      )
      
      ## DIAMOND ------------------------------------------------------
    } else if (symbol == "diamond") {
      
      graphics::polygon(
        c(x[i], x[i] + rx, x[i], x[i] - rx),
        c(y[i] + ry, y[i], y[i] - ry, y[i]),
        col = NA,
        border = col,
        lwd = lwd
      )
      
      ## TRIANGLE -----------------------------------------------------
    } else if (symbol == "triangle") {
      
      graphics::polygon(
        c(x[i], x[i] + rx, x[i] - rx),
        c(y[i] + ry, y[i] - ry, y[i] - ry),
        col = NA,
        border = col,
        lwd = lwd
      )
      
      } else {
      stop("Supported symbols are 'rectangle', 'ellipse', 'diamond', 'triangle', 'cross', 'plus', 'star'")
    }
  }
  invisible(NULL)
}