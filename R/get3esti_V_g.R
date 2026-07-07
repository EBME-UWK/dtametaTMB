#' @keywords internal
#' @noRd
get3esti_V_g <- function(beta_fix, theta, V_full) {
  ## transformed parameters
  sigma_A <- exp(theta[1])
  sigma_B <- exp(theta[2])
  rho     <- theta[3] / sqrt(1 + theta[3]^2)
  ## 
  g <- c(beta_fix,
         exp(2 * theta[1]),
         exp(2 * theta[2]),
         rho * sigma_A * sigma_B)

  
  ## --- Jacobian J = dg / d(beta_fix, theta) ---
  ## parameter order assumed in V_full:
  ## (beta_fix[1], beta_fix[2], theta[1], theta[2], theta[3])
  J <- matrix(0, nrow = 5, ncol = 5)
  
  ## mu_A, mu_B
  J[1, 1] <- 1
  J[2, 2] <- 1
  
  ## sigma2_A = exp(2*theta1)
  J[3, 3] <- 2 * exp(2 * theta[1])
  
  ## sigma2_B = exp(2*theta2)
  J[4, 4] <- 2 * exp(2 * theta[2])
  
  ## sigma_AB = rho * sigma_A * sigma_B
  sigma_AB <- rho * sigma_A * sigma_B
  
  ## d sigma_AB / d theta1 = sigma_AB
  J[5, 3] <- sigma_AB
  
  ## d sigma_AB / d theta2 = sigma_AB
  J[5, 4] <- sigma_AB
  
  ## d sigma_AB / d theta3
  ## rho = theta3 / sqrt(1 + theta3^2)
  ## drho/dtheta3 = 1 / (1 + theta3^2)^(3/2)
  J[5, 5] <- sigma_A * sigma_B / (1 + theta[3]^2)^(3/2)
  
  ## --- delta-method covariance ---
  V_g <- J %*% V_full %*% t(J)
  colnames(V_g) <- rownames(V_g) <- names(g)
  ## --- summary table ---
  esti <- data.frame(
    Estimate  = g,
    Std_Error = sqrt(diag(V_g))
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
