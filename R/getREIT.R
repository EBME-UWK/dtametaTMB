#' @keywords internal
#' @noRd
getREIT<- function(Lambda,
                   Theta,
                   beta,
                   sigma2_alpha,
                   sigma2_theta) {
  b <- exp(beta / 2)
  
  reit <- data.frame(
    mu_A.sens     = b^(-1) * (Theta + 0.5 * Lambda),
    mu_B.spec     = -b * (Theta - 0.5 * Lambda),
    sigma2_A.sens = b^(-2) * (sigma2_theta + 0.25 * sigma2_alpha),
    sigma2_B.spec = b^( 2) * (sigma2_theta + 0.25 * sigma2_alpha),
    sigma_AB      = -(sigma2_theta - 0.25 * sigma2_alpha),
    row.names     = "Estimate (recovered)"
  )
  return(reit)
}