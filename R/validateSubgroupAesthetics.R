#' @keywords internal
#' @noRd
validateSubgroupAesthetics <- function(col,symbol,nsub) {
  
  available_symbols <- c("rectangle","ellipse","diamond","triangle","cross","plus","star")
  
  # colour
  if (is.null(col)) {col <- grDevices::rainbow(n = nsub)} else {
    if (!is.character(col) || !length(col) %in% c(1L, nsub)) {
      stop("'col' must be NULL, a single colour, or a character ",
           "vector with one colour for each subgroup (expected length: ",
           nsub,").")}
    if (length(col) == 1L) { col <- rep(col, nsub) }
    invalid_col <- vapply(col,function(z) { tryCatch(  { grDevices::col2rgb(z); FALSE }, error = function(e) TRUE  ) }, logical(1) )
    if (any(invalid_col)) {
      stop("Invalid colour value(s) in 'col': ",
           paste(shQuote(col[invalid_col]),collapse = ", "),".")}
  }
  # symbol
  if (is.null(symbol)) {symbol <- rep(available_symbols,length.out = nsub)} else {
    if (!is.character(symbol) || !length(symbol) %in% c(1L, nsub)) {
      stop("'symbol' must be NULL, a single symbol, or a character ",
           "vector with one symbol for each subgroup (expected length: ",
           nsub,").")}
    if (length(symbol) == 1L) { symbol <- rep(symbol, nsub) }
    invalid_symbol <- setdiff(unique(symbol),available_symbols)
    if (length(invalid_symbol) > 0L) {
      stop("Unknown symbol value(s): ",
           paste(shQuote(invalid_symbol),collapse = ", "),
           ". Available symbols are: ",
           paste(shQuote(available_symbols),collapse = ", "),".")}
  }
  
  list(col = col,
       col2 = grDevices::adjustcolor(col,alpha.f = 0.6),
       symbols2 = symbol)
}