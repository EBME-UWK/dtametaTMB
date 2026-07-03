#' @keywords internal
#' @noRd
getLRDOR <- function(lsens, lspec, S, conflevel) {
  sens <- plogis(lsens)
  spec <- plogis(lspec)
  ## point estimates
  DOR <- exp(lsens + lspec)
  LRp <- sens / (1 - spec)
  LRn <- (1 - sens) / spec
  ## delta-method SEs on log scale
  # log(DOR) = lsens + lspec
  gDOR <- c(1, 1)
  se.logDOR <- sqrt(drop(t(gDOR) %*% S %*% gDOR))
  # log(LR+) = log(sens) - log(1 - spec)
  gLRp <- c(1 / (1 + exp(lsens)),
            exp(lspec) / (1 + exp(lspec)))
  se.logLRp <- sqrt(drop(t(gLRp) %*% S %*% gLRp))
  # log(LR-) = log(1 - sens) - log(spec)
  gLRn <- c(-exp(lsens) / (1 + exp(lsens)),
            -1 / (1 + exp(lspec)))
  se.logLRn <- sqrt(drop(t(gLRn) %*% S %*% gLRn))
  ## output
  qq    <- stats::qnorm(1-(1-conflevel)/2)
  lrdor <- data.frame(
    Estimate  = c(DOR, LRp, LRn),
    conflevel = conflevel,
    CI_Lower  = c(exp(log(DOR) - qq * se.logDOR),
                  exp(log(LRp) - qq * se.logLRp),
                  exp(log(LRn) - qq * se.logLRn)),
    CI_Upper  = c(exp(log(DOR) + qq * se.logDOR),
                  exp(log(LRp) + qq * se.logLRp),
                  exp(log(LRn) + qq * se.logLRn)),
    row.names = c("DOR", "LR+", "LR-")
  )
  return(lrdor)
}