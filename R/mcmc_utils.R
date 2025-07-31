comp_lgamma_lk_R <- function(gamma_lk, FF, Gamma, Z_sum, delta, l, k){
  Gamma_lk <- Gamma
  Gamma_lk[l,k] <- gamma_lk
  
  FF_Gamma_lk <- t(t(FF%*%Gamma_lk) + delta)
  FFG_log_exp_sums <- log(rowSums(exp(FF_Gamma_lk)))
  
  lgamma_lk <- sum(Z_sum*(FF_Gamma_lk - FFG_log_exp_sums)) - 1/2*gamma_lk^2
  
  lgamma_lk
}

sample_gamma_lk_R <- function(FF, Gamma, Z_sum, delta, gamma_sd, l, k){
  gamma_lk <- Gamma[l,k]
  Gamma_star <- Gamma
  gamma_lk_star <- rnorm(1, gamma_lk, gamma_sd)
  Gamma_star[l,k] <- gamma_lk_star
  
  log_gamma_lk <- comp_lgamma_lk_R(gamma_lk, FF, Gamma, Z_sum, delta, l, k)
  log_gamma_lk_star <- comp_lgamma_lk_R(gamma_lk_star, FF, Gamma, Z_sum, delta, l, k)
  
  log_u <- log(runif(1))
  log_alpha <- log_gamma_lk_star - log_gamma_lk

  if (log_alpha > log_u){
    gamma_lk <- gamma_lk_star
  }
  
  gamma_lk
}

sample_Gamma_R <- function(FF, Gamma, Z_sum, delta, gamma_sd){
  L <- nrow(Gamma)
  K <- ncol(Gamma)
  
  for (l in 1:L){
    for (k in 1:(K-1)){
      Gamma[l,k] <- sample_gamma_lk_R(FF, Gamma, Z_sum, delta, gamma_sd, l, k)
    }
  }
  
  Gamma
}

comp_lf_d_R <- function(f_d, Y_d, Z_sum_d, Phi, Gamma, R_inv, mu, sigma2, delta){
  f_d <- rbind(f_d)
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

sample_f_d_R <- function(f_d, Y_d, Z_sum_d, Phi, Gamma, R_inv, 
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
    FF[d,] <- sample_f_d_R(f_d, Y[d,], Z_sum[d,], Phi, Gamma, R_inv, 
                           mu, sigma2, delta, f_Sigma)
  }
  
  FF
}

comp_ldelta_k_R <- function(delta_k, FF, Gamma, Z_sum, delta, k){
  delta_k_vec <- delta
  delta_k_vec[k] <- delta_k
  
  FF_Gamma_k <- t(t(FF%*%Gamma) + delta_k_vec)
  FFG_log_exp_sums <- log(rowSums(exp(FF_Gamma_k)))
  
  ldelta_k <- sum(Z_sum*(FF_Gamma_k - FFG_log_exp_sums)) - 1/2*delta_k^2
  
  ldelta_k
}

sample_delta_k <- function(FF, Gamma, Z_sum, delta, delta_sd, k){
  delta_k <- delta[k]
  delta_k_star <- rnorm(1, delta_k, delta_sd)
  
  log_delta_k <- comp_ldelta_k_R(delta_k, FF, Gamma, Z_sum, delta, k)
  log_delta_k_star <- comp_ldelta_k_R(delta_k_star, FF, Gamma, Z_sum, delta, k)
  
  log_u <- log(runif(1))
  log_alpha <- log_delta_k_star - log_delta_k
  
  if (log_alpha > log_u){
    delta_k <- delta_k_star
  }
  
  delta_k
}

sample_delta <- function(FF, Gamma, Z_sum, delta, delta_sd){
  K <- length(delta)
  
  for (k in 1:(K-1)){
    delta[k] <- sample_delta_k(FF, Gamma, Z_sum, delta, delta_sd, k)
  }
  
  delta
}

comp_Z_post_R <- function(W, FF, Gamma, Psi, delta){
  D <- length(W)
  Z_post <- vector("list", D)
  FF_Gamma <- t(t(FF%*%Gamma) + delta)
  Theta <- exp(FF_Gamma)/rowSums(exp(FF_Gamma))
  
  for (d in 1:D){
    Psi_d <- Psi[W[[d]],]
    Z_post_unnorm_d <- t(Psi_d)*Theta[d,]
    Z_post_d <- t(Z_post_unnorm_d)/colSums(Z_post_unnorm_d)
    Z_post[[d]] <- Z_post_d
  }
  
  Z_post
}

sample_Z_R <- function(W, FF, Gamma, Psi, delta){
  D <- length(W)
  K <- ncol(Gamma)
  Z <- vector("list", D)
  
  Z_post <- comp_Z_post_R(W, FF, Gamma, Psi, delta)

  for (d in 1:D){
    Z[[d]] <- apply(Z_post[[1]], 1, function(p){
      sample(1:K, 1, FALSE, p)
    })
  }  
  
  Z
}
