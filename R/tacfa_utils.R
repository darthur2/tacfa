sum_Z_R <- function(Z){
  D <- length(Z)
  K <- ncol(Z[[1]])
  Z_sum <- matrix(0, D, K)
  
  for (d in 1:D){
    Z_sum[d,] <- colSums(Z[[d]])
  }
  
  Z_sum
}

comp_lgamma_lk_R <- function(gamma_lk, FF, Gamma, Z_sum, delta, l, k){
  Gamma_lk <- Gamma
  Gamma_lk[l,k] <- gamma_lk
  
  FF_Gamma_lk <- t(t(FF%*%Gamma_lk) + delta)
  FFG_log_exp_sums <- log(rowSums(exp(FF_Gamma_lk)))
  
  lgamma_lk <- sum(Z_sum*(FF_Gamma_lk - FFG_log_exp_sums)) - 1/2*gamma_lk^2

  lgamma_lk
}

comp_lf_d_R <- function(f_d, Y_d, Z_sum_d, Phi, Gamma, R_inv, mu, sigma2, delta){
  Phi_sigma <- t(Phi)*1/sqrt(sigma2)
  f_d_Sigma_inv <- t(Phi_sigma)%*%Phi_sigma + R_inv
  f_d_Mu <- colSums(t(Phi)*((Y_d-mu)/sigma2))
  
  f_d_Gamma <- f_d%*%Gamma + delta
  theta_d <- exp(f_d_Gamma)/sum(exp(f_d_Gamma))
  Z_sum_log_theta_d <- sum(Z_sum_d*log(theta_d))
  
  lf_d <- as.numeric(-1/2*(f_d%*%f_d_Sigma_inv%*%t(f_d) -
                             2*f_d%*%f_d_Mu) +
    Z_sum_log_theta_d)
  
  lf_d
}

