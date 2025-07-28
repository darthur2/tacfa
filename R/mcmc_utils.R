sample_gamma_lk <- function(FF, Gamma, Z_sum, delta, gamma_sd, l, k){
  gamma_lk <- Gamma[l,k]
  Gamma_star <- Gamma
  gamma_lk_star <- rnorm(1, gamma_lk, gamma_sd)
  Gamma_star[l,k] <- gamma_lk_star
  
  log_gamma_lk <- comp_lgamma_lk_cpp(gamma_lk, FF, Gamma, Z_sum, delta, l, k)
  log_gamma_lk_star <- comp_lgamma_lk_cpp(gamma_lk_star, FF, Gamma, Z_sum, delta, l, k)
  
  log_u <- log(runif(1))
  log_alpha <- log_gamma_lk_star - log_gamma_lk

  if (log_alpha > log_u){
    gamma_lk <- gamma_lk_star
  }
  
  gamma_lk
}

sample_Gamma <- function(FF, Gamma, Z_sum, delta, gamma_sd){
  L <- nrow(Gamma)
  K <- ncol(Gamma)
  
  for (l in 1:L){
    for (k in 1:(K-1)){
      Gamma[l,k] <- sample_gamma_lk(FF, Gamma, Z_sum, delta, gamma_sd, l, k)
    }
  }
  
  Gamma
}

sample_f_d <- function(f_d, Y_d, Z_sum_d, Phi, Gamma, R_inv, 
                       mu, sigma2, delta, f_Sigma){
  f_d_star <- rbind(mvrnorm(1, f_d, f_Sigma))
  
  log_f_d <- comp_lf_d_R(f_d, Y_d, Z_sum_d, Phi, Gamma, R_inv, mu, sigma2, delta)
  log_f_d_star <- comp_lf_d_R(f_d_star, Y_d, Z_sum_d, Phi, Gamma, R_inv, mu, sigma2, delta)
  
  log_u <- log(runif(1))
  log_alpha <- log_f_d_star - log_f_d
  
  if (log_alpha > log_u){
    f_d <- f_d_star
  }
  
  f_d
}

sample_FF <- function(FF, Y, Z_sum, Phi, Gamma, R_inv, mu, sigma2, delta, f_Sigma){
  D <- nrow(FF)
  
  for (d in 1:D){
    f_d <- rbind(FF[d,])
    FF[d,] <- sample_f_d(f_d, Y[d,], Z_sum[d,], Phi, Gamma, R_inv, 
                         mu, sigma2, delta, f_Sigma)
  }
  
  FF
}
