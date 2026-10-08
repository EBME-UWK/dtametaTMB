#' Draw Study-Level Confidence Intervals on an SROC Plot
#'
#' @keywords internal
#' @importFrom graphics segments
#' @noRd
drawStudyPointCI <- function(
    x,
    y,
    sens_lower,
    sens_upper,
    spec_lower,
    spec_upper,
    symbol,
    scale,
    cex.x,
    cex.y,
    col = "grey70",
    lwd = 1,
    lty = 1) {
  
  n <- length(x)
  
  required_lengths <- c(
    length(y),
    length(sens_lower),
    length(sens_upper),
    length(spec_lower),
    length(spec_upper),
    length(cex.x),
    length(cex.y)
  )
  
  if (any(required_lengths != n)) {
    stop("All study-level graphical inputs must have the same length.")
  }
  
  if (length(symbol) == 1L) {
    symbol <- rep(symbol, n)
  }
  
  if (length(symbol) != n) {
    stop("'symbol' must have length one or one value per study.")
  }
  
  ext <- getStudySymbolExtent(symbol = symbol,
                              scale = scale,
                              cex.x = cex.x,
                              cex.y = cex.y)
  
  # Convert specificity limits to false-positive-rate coordinates.
  x_lower <- 1 - spec_upper
  x_upper <- 1 - spec_lower
  
  valid <- is.finite(x) &
    is.finite(y) &
    is.finite(x_lower) &
    is.finite(x_upper) &
    is.finite(sens_lower) &
    is.finite(sens_upper)
  
  for (i in which(valid)) {
      # Horizontal confidence interval, left side of symbol.
      left_end <- x[i] - ext$rx[i]
      
      if (x_lower[i] < left_end) {
        graphics::segments(
          x0 = x_lower[i],
          y0 = y[i],
          x1 = left_end,
          y1 = y[i],
          col = col,
          lwd = lwd,
          lty = lty
        )
      }
      
      # Horizontal confidence interval, right side of symbol.
      right_start <- x[i] + ext$rx[i]
      
      if (right_start < x_upper[i]) {
        graphics::segments(
          x0 = right_start,
          y0 = y[i],
          x1 = x_upper[i],
          y1 = y[i],
          col = col,
          lwd = lwd,
          lty = lty
        )
      }
      
      # Vertical confidence interval, below symbol.
      lower_end <- y[i] - ext$ry[i]
      
      if (sens_lower[i] < lower_end) {
        graphics::segments(
          x0 = x[i],
          y0 = sens_lower[i],
          x1 = x[i],
          y1 = lower_end,
          col = col,
          lwd = lwd,
          lty = lty
        )
      }
      
      # Vertical confidence interval, above symbol.
      upper_start <- y[i] + ext$ry[i]
      
      if (upper_start < sens_upper[i]) {
        graphics::segments(
          x0 = x[i],
          y0 = upper_start,
          x1 = x[i],
          y1 = sens_upper[i],
          col = col,
          lwd = lwd,
          lty = lty
        )
      }
      
    }
  invisible(NULL)
}


#' Calculate Half-Extents of SROC Study Symbols
#'
#' @keywords internal
#' @noRd
getStudySymbolExtent <- function(symbol, scale, cex.x, cex.y) {
  
  if (length(symbol) == 1L) {
    symbol <- rep(symbol, length(cex.x))
  }
  
  if (length(symbol) != length(cex.x) ||
      length(cex.x) != length(cex.y)) {
    stop(
      "'symbol', 'cex.x', and 'cex.y' must have compatible lengths."
    )
  }
  
  shape_factor <- rep(1, length(symbol))
  
  shape_factor[symbol == "ellipse"] <- 2 / sqrt(pi)
  
  shape_factor[symbol %in% c("diamond", "triangle")] <- sqrt(2)
  
  data.frame(
    rx = scale * cex.x * shape_factor,
    ry = scale * cex.y * shape_factor
  )
}