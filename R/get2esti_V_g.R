#' @keywords internal
#' @noRd
get2esti_V_g <- function(beta_fix, theta, V_full) {
  ## transformed parameters
  sigma_A <- exp(theta[1])
  sigma_B <- exp(theta[2])
  rho     <- theta[3] / sqrt(1 + theta[3]^2)
  ## 
  g <- c(beta_fix,
    sigma2_A.sens = exp(2 * theta[1]),
    sigma2_B.spec = exp(2 * theta[2]),
    sigma_AB      = rho * sigma_A * sigma_B
  )
 
  ## --- Jacobian J = dg / d(beta_fix, theta) ---
  ## parameter order assumed in V_full:
  ## 
  lb <- length(beta_fix)
  lg <- length(g)
  J  <- matrix(0, nrow = lg, ncol = lg)
  
  ## mu_A, mu_B, vu
  for(i in seq_len(lb)) { J[i,i] <- 1 }
  
  ## sigma2_A = exp(2*theta1)
  J[lb+1, lb+1] <- 2 * exp(2 * theta[1])
  
  ## sigma2_B = exp(2*theta2)
  J[lb+2, lb+2] <- 2 * exp(2 * theta[2])
  
  ## sigma_AB = rho * sigma_A * sigma_B
  sigma_AB <- rho * sigma_A * sigma_B
  
  ## d sigma_AB / d theta1 = sigma_AB
  J[lb+3, lb+1] <- sigma_AB
  
  ## d sigma_AB / d theta2 = sigma_AB
  J[lb+3, lb+2] <- sigma_AB
  
  ## d sigma_AB / d theta3
  ## rho = theta3 / sqrt(1 + theta3^2)
  ## drho/dtheta3 = 1 / (1 + theta3^2)^(3/2)
  J[lb+3, lb+3] <- sigma_A * sigma_B / (1 + theta[3]^2)^(3/2)
  colnames(J) <- rownames(J) <- names(g)
  ## --- delta-method covariance ---
  V_g <- J %*% V_full %*% t(J)
  
  ## --- summary table ---
  esti <- data.frame(
    Estimate  = g,
    Std_Error = sqrt(diag(V_g)),
    row.names = names(g)
  )
  #esti$CI_Lower <- with(esti,Estimate-qq*Std_Error)
  #esti$CI_Upper <- with(esti,Estimate+qq*Std_Error)
  
  ## return both
  list(
    esti = esti,
    V_g  = V_g,
    J    = J
  )
}