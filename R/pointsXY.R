#' Draw Study Symbols with Separate Horizontal and Vertical Dimensions
#'
#' Internal graphical helper for drawing study symbols centred at supplied
#' plotting coordinates. The \code{cex.x} and \code{cex.y} arguments specify
#' horizontal and vertical dimension measures, which are multiplied by
#' \code{scale} to obtain the reference rectangle half-width and half-height.
#'
#' For \code{symbol = "rectangle"}, the reference half-width and half-height
#' are used directly. Ellipses, diamonds, and triangles are rescaled so that
#' their enclosed area equals the area of the corresponding reference
#' rectangle. The ellipse semiaxes are multiplied by
#' \eqn{2 / \sqrt{\pi}}, whereas the horizontal and vertical half-extents
#' of diamonds and triangles are multiplied by \eqn{\sqrt{2}}.
#'
#' The equal-area adjustment preserves the ratio between the horizontal and
#' vertical dimensions. Consequently, horizontal size continues to represent
#' specificity-side information and vertical size continues to represent
#' sensitivity-side information.
#'
#' The line-based symbols \code{"plus"}, \code{"cross"}, and \code{"star"}
#' do not enclose a two-dimensional area and therefore use the unadjusted
#' reference half-width and half-height.
#'
#' @param x Numeric vector of horizontal plotting coordinates.
#' @param y Numeric vector of vertical plotting coordinates.
#' @param symbol Character string specifying the symbol. Supported values are
#'   \code{"rectangle"}, \code{"ellipse"}, \code{"diamond"},
#'   \code{"triangle"}, \code{"cross"}, \code{"plus"}, and \code{"star"}.
#' @param scale Positive numeric multiplier controlling the overall symbol
#'   size.
#' @param cex.x Numeric vector of horizontal dimension measures. Must have
#'   the same length as \code{x} and \code{y}.
#' @param cex.y Numeric vector of vertical dimension measures. Must have
#'   the same length as \code{x} and \code{y}.
#' @param col Symbol border or line colour.
#' @param lwd Symbol border or line width.
#'
#' @return Invisibly returns \code{NULL}. The function is called for its
#'   graphical side effect.
#'
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
            length(x) == length(cex.y) &
            length(scale) == 1L &
            is.numeric(scale) &
            is.finite(scale) &
            scale > 0)
  
  for (i in seq_along(x)) {
    ## rx and ry define the half-width and half-height of the
    ## reference rectangle. Closed non-rectangular symbols are
    ## rescaled below to preserve the same enclosed area.
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
      
      theta <- seq(0, 2*pi, length.out = 200)
      
      ## Preserve the area of the corresponding reference rectangle.
      shape_factor <- 2 / sqrt(pi)
      
      graphics::polygon(
        x[i] + shape_factor * rx * cos(theta),
        y[i] + shape_factor * ry * sin(theta),
        col = NA,
        border = col,
        lwd = lwd
      )
      
      ## RECTANGLE ---------------------------------------------------
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
      
      ## Preserve the area of the corresponding reference rectangle.
      shape_factor <- sqrt(2)
      
      graphics::polygon(
        c(x[i], 
          x[i] + shape_factor * rx, 
          x[i], 
          x[i] - shape_factor * rx),
        c(y[i] + shape_factor * ry, 
          y[i], 
          y[i] - shape_factor * ry, 
          y[i]),
        col = NA,
        border = col,
        lwd = lwd
      )
      
      ## TRIANGLE -----------------------------------------------------
    } else if (symbol == "triangle") {
      
      ## Preserve the area of the corresponding reference rectangle.
      shape_factor <- sqrt(2)
      
      graphics::polygon(
        c(x[i], 
          x[i] + shape_factor * rx, 
          x[i] - shape_factor * rx),
        c(y[i] + shape_factor * ry, 
          y[i] - shape_factor * ry, 
          y[i] - shape_factor * ry),
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