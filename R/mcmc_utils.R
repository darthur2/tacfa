sum_Z_R <- function(Z, K) {
  D <- length(Z)
  Z_sum <- matrix(0, nrow = D, ncol = K)
  
  for (d in seq_len(D)) {
    z_d <- Z[[d]]
    for (k in z_d) {
      Z_sum[d,k] <- Z_sum[d,k] + 1
    }
  }
  
  Z_sum
}

rdirichlet_R <- function(n, shape){
  m <- length(shape)
  draws <- matrix(0, m, n)
  
  for (i in 1:n){
    gamma_draws_i <- rgamma(m, shape)
    draws[,i] <- gamma_draws_i/sum(gamma_draws_i)
  }
  
  draws
}

comp_WZ_counts_R <- function(W, Z, V, K) {
  WZ_counts <- matrix(0, nrow = V, ncol = K)
  
  for (d in seq_along(W)) {
    w_d <- W[[d]]
    z_d <- Z[[d]]
    for (i in seq_along(w_d)) {
      word <- w_d[i]
      topic <- z_d[i]
      WZ_counts[word, topic] <- WZ_counts[word, topic] + 1
    }
  }
  
  WZ_counts  
}

rmvnorm_R <- function(n, Mu, Sigma){
  Q <- chol(Sigma)
  L <- length(Mu)
  
  X <- matrix(rnorm(n*L), n, L)
  Y <- t(t(X%*%Q) + Mu)
  
  Y
}

rtnorm_R <- function(n, mu, sigma, lower, upper){
  p_lower <- pnorm(lower, mu, sigma)
  p_upper <- pnorm(upper, mu, sigma)
  
  u <- runif(n, p_lower, p_upper)
  draws <- qnorm(u)*sigma + mu
  
  draws
}

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
    for (k in 1:(K)){
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
  f_d_star <- rbind(rmvnorm_R(1, f_d, f_Sigma))
  
  log_f_d <- comp_lf_d_R(f_d, Y_d, Z_sum_d, Phi, Gamma, R_inv, mu, sigma2, delta)
  log_f_d_star <- comp_lf_d_R(f_d_star, Y_d, Z_sum_d, Phi, Gamma, R_inv, mu, sigma2, delta)
  
  log_u <- log(runif(1))
  log_alpha <- log_f_d_star - log_f_d
  
  if (log_alpha > log_u){
    f_d <- f_d_star
  }
  
  f_d
}

sample_FF_R <- function(FF, Y, Z_sum, Phi, Gamma, R_inv, mu, sigma2, delta, f_Sigma){
  D <- nrow(FF)
  
  for (d in 1:D){
    f_d <- FF[d,]
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

sample_delta_k_R <- function(FF, Gamma, Z_sum, delta, delta_sd, k){
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

sample_delta_R <- function(FF, Gamma, Z_sum, delta, delta_sd){
  K <- length(delta)
  
  for (k in 1:(K-1)){
    delta[k] <- sample_delta_k_R(FF, Gamma, Z_sum, delta, delta_sd, k)
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
    Z[[d]] <- apply(Z_post[[d]], 1, function(p){
      sample(1:K, 1, FALSE, p)
    })
  }  
  
  Z
}

sample_phi_lj_R <- function(Y, FF, mu, sigma2, Phi, l, j){
  FF_Phi_lj <- cbind(FF[,-l])%*%rbind(Phi[-l,j])
  
  phi_lj_var <- 1/(1/sigma2[j]*sum(FF[,l]^2) + 1)
  phi_lj_mu <- 1/sigma2[j]*sum((Y[,j] - mu[j] - FF_Phi_lj)*FF[,l])

  phi_lj <- rtnorm_R(1, phi_lj_mu*phi_lj_var, sqrt(phi_lj_var), 0, Inf)
  
  phi_lj
}

sample_Phi_R <- function(Y, FF, mu, sigma2, Phi, Phi_indices){
  
  for (l in 1:L){
    for (j in Phi_indices[[l]]){
      Phi[l,j] <- sample_phi_lj_R(Y, FF, mu, sigma2, Phi, l, j)
    }
  }
  
  Phi
}

sample_sigma2_j_R <- function(Y, mu, FF, Phi, a_sigma2, b_sigma2, j){
  D <- nrow(FF)

  FF_Phi_j <- (FF%*%Phi)[,j]

  a_sigma2_star <- D/2 + a_sigma2
  b_sigma2_star <- 1/2*sum((Y[,j] - mu[j] - FF_Phi_j)^2) + b_sigma2

  sigma2_j <- 1/rgamma(1, a_sigma2_star, b_sigma2_star)

  sigma2_j
}

sample_sigma2_R <- function(Y, mu, FF, Phi, a_sigma2, b_sigma2){
  
  J <- ncol(Phi)
  sigma2 <- rep(0, J)
  
  for (j in 1:J){
    sigma2[j] <- sample_sigma2_j_R(Y, mu, FF, Phi, a_sigma2, b_sigma2, j)
  }
  
  sigma2
}

sample_mu_j_R <- function(Y, sigma2, FF, Phi, sigma2_mu, j){
  D <- nrow(FF)
  FF_Phi_j <- (FF%*%Phi)[,j]
  
  mu_j_var <- 1/(D/sigma2[j] + 1/sigma2_mu)
  mu_j_mu <- 1/sigma2[j]*sum(Y[,j] - FF_Phi_j)
  
  mu_j <- rnorm(1, mu_j_mu*mu_j_var, sqrt(mu_j_var))
  
  mu_j
}

sample_mu_R <- function(Y, sigma2, FF, Phi, sigma2_mu){
  J <- ncol(Phi)
  mu <- rep(0, J)
    
  for (j in 1:J){
    mu[j] <- sample_mu_j_R(Y, sigma2, FF, Phi, sigma2_mu, j)
  }
  
  mu
}

sample_Psi_R <- function(W, Z, beta){
  V <- nrow(beta)
  K <- ncol(beta)
  
  Psi <- matrix(0, V, K)
  WZ_counts <- comp_WZ_counts_R(W, Z, V, K)
  beta_star <- WZ_counts + beta
  
  for (k in 1:K){
    Psi[,k] <- rdirichlet_R(1, beta_star[,k])
  }
  
  Psi
}
