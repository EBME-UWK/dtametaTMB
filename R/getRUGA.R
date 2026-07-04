#' @keywords internal
#' @noRd
getRUGA <- function(lsens,
                    lspec,
                    sigma_a,
                    sigma_b,
                    sigma_ab) {
  ruga <- data.frame(
    Lambda   = (((sigma_b/sigma_a)**0.5) * lsens) + ((sigma_a/sigma_b)**0.5 *lspec),
    Theta    = 0.5*((((sigma_b/sigma_a)**0.5 )*lsens) - (((sigma_a/sigma_b)**0.5) *lspec)),
    beta     = log(sigma_b/sigma_a),
    sigma2_alpha = 2*((sigma_a*sigma_b) + sigma_ab),
    sigma2_theta = 0.5*((sigma_a*sigma_b) - sigma_ab),
    row.names="Estimate (recovered)")
  return(ruga)
}