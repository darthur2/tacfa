tacfaSM <- function(Y,
                    W,
                    V,
                    phi_factor,
                    phi_fixed,
                    gamma_factor,
                    gamma_fixed,
                    gamma_fixed_value,
                    topic_baseline,
                    anchor_words,
                    anchor_prob,
                    nonanchor_prob,
                    mu_unit,
                    mu_0 = 0,
                    sigma_mu2 = 1,
                    mu_phi = 0,
                    sigma_phi2 = 1,
                    a_sigma2 = 2,
                    b_sigma2 = 2,
                    a_rho = 1,
                    b_rho = 1,
                    mu_gamma = 0,
                    sigma_gamma2 = 1,
                    mu_tau = 0,
                    sigma_tau2 = 1,
                    beta_const = 0.01,
                    delta = 1e-05,
                    target_accept = 0.574,
                    kappa = 0.6,
                    niter = 5000,
                    nburn = 1000,
                    warmup = 250){
  N <- nrow(Y)
  J <- ncol(Y)
  L <- max(phi_factor)
  K <- length(gamma_factor) + L + 1
  
  mu <- sapply(1:J, function(j) rtnorm_lower_single(mu_phi, sigma_phi2))
  
  phi <- rnorm(J-L, mu_phi, sigma_phi2)
  Phi <- make_Phi(phi, phi_factor, phi_fixed)
  
  sigma2 <- 1/rgamma(J, a_sigma2, b_sigma2)
  
  rho <- rbeta(1, a_rho, b_rho)
  log_step_rho <- log(0.1)
  R <- make_R(rho, L)
  chol_R <- chol(R)
  R_inv <- chol2inv(chol_R)
  
  FF <- matrix(rnorm(N*L), N, L)
  log_step_FF <- log(0.1)
  G_FF <- R
  chol_G_FF <- chol(G_FF)
  G_inv_FF <- chol2inv(chol_G_FF)
  FF_mean <- colMeans(FF)
  C_FF <- matrix(0, L, L)
  
  gamma <- rnorm(K-L-1, mu_gamma, sigma_gamma2)
  gamma_args <- list(gamma = gamma,
                     gamma_factor = gamma_factor,
                     gamma_fixed = gamma_fixed,
                     gamma_fixed_value = gamma_fixed_value,
                     topic_baseline = topic_baseline)
  Gamma <- do.call(make_Gamma, gamma_args)
  log_step_gamma <- log(0.1)
  G_inv_gamma <- G_gamma <- diag(K-L-1)
  gamma_mean <- rep(0, K-L-1)
  C_gamma <- matrix(0, K-L-1, K-L-1)
  
  tau <- rnorm(K-1, mu_tau, sigma_tau2)
  tau_args <- list(tau = tau,
                   topic_baseline = topic_baseline)
  tau_full <- make_tau_full(tau, topic_baseline)
  log_step_tau <- log(0.1)
  G_inv_tau <- G_tau <- diag(K-1)
  tau_mean <- rep(0, K-1)
  C_tau <- matrix(0, K-1, K-1)
  
  Theta <- comp_Theta(tau_full, FF, Gamma)
  
  M <- sapply(W, length)
  Z <- make_Z(Theta, M)
  n_ik <- comp_n_ik(Z, K)
  
  n_vk <- comp_n_vk(W, Z, V, K)
  
  beta <- rep(beta_const, V-K)
  Psi_full <- matrix(0, V, K)
  Psi_full[cbind(anchor_words, 1:K)] <- anchor_prob - nonanchor_prob
  Psi_full[anchor_words, 1:K] <- Psi_full[anchor_words, 1:K] + nonanchor_prob
  
  for (k in 1:K){
    dir_sample <- rdirichlet_single(beta)
    Psi_full[-anchor_words,k] <- (1-anchor_prob-(K-1)*nonanchor_prob)*dir_sample
  }
  
  mu_save <- matrix(0, niter-nburn, J)
  phi_save <- matrix(0, niter-nburn, J-L)
  sigma2_save <- matrix(0, niter-nburn, J)
  FF_save <- matrix(0, niter-nburn, N*L)
  gamma_save <- matrix(0, niter-nburn, K-L-1)
  tau_save <- matrix(0, niter-nburn, K-1)

  for (iter in 1:niter){
    mu <- sample_mu(Y, FF, Phi, sigma2, mu_0, sigma_mu2)
    
    phi <- sample_phi(Y, FF, mu, sigma2, mu_phi, sigma_phi2, phi_fixed, phi_factor)
    Phi <- make_Phi(phi, phi_factor, phi_fixed)
    
    sigma2 <- sample_sigma2(Y, FF, Phi, mu, a_sigma2, b_sigma2)
    
    FF_out <- sample_FF(FF, Y, n_ik, mu, R_inv, sigma2, Phi, Gamma, tau_full,
                        log_step_FF, iter, warmup, G_FF, G_inv_FF, FF_mean,
                        C_FF, delta, target_accept, kappa)
    FF <- FF_out$FF
    log_step_FF <- FF_out$log_step
    G_FF <- FF_out$G_FF
    G_inv_FF <- FF_out$G_inv_FF
    FF_mean <- FF_out$FF_mean
    C_FF <- FF_out$C_FF
    
    rho_out <- sample_rho(rho, R, FF, a_rho, b_rho, log_step_rho, iter, nburn,
                          target_accept, kappa)
    rho <- rho_out$rho
    R <- rho_out$R
    chol_R <- chol(R)
    R_inv <- chol2inv(chol_R)
    
    gamma_out <- sample_gamma(FF, gamma, n_ik, tau_full, mu_gamma, sigma_gamma2,
                              gamma_args, log_step_gamma, iter, warmup, G_gamma,
                              G_inv_gamma, gamma_mean, C_gamma, delta,
                              target_accept, kappa)
    gamma <- gamma_out$gamma
    Gamma <- gamma_out$Gamma
    log_step_gamma <- gamma_out$log_step
    G_gamma <- gamma_out$G_gamma
    G_inv_gamma <- gamma_out$G_inv_gamma
    gamma_mean <- gamma_out$gamma_mean
    C_gamma <- gamma_out$C_gamma
    
    FF_Gamma <- FF%*%Gamma
    
    tau_out <- sample_tau(FF_Gamma, tau, n_ik, mu_tau, sigma_tau2, tau_args,
                          log_step_tau, iter, warmup, G_tau, G_inv_tau, tau_mean,
                          C_tau, delta, target_accept, kappa)
    tau <- tau_out$tau
    tau_full <- tau_out$tau_full
    log_step_tau <- tau_out$log_step
    G_tau <- tau_out$G_tau
    G_inv_tau <- tau_out$G_inv_tau
    tau_mean <- tau_out$tau_mean
    C_tau <- tau_out$C_tau
    
    Z <- sample_Z(W, FF, tau_full, Gamma, Psi_full)
    n_ik <- comp_n_ik(Z, K)
    n_vk <- comp_n_vk(W, Z, V, K)
    
    Psi <- sample_Psi(beta, n_vk, anchor_words, anchor_prob, nonanchor_prob)
    Psi_full <- make_Psi_full(Psi_full, Psi, anchor_words)
    
    if (iter > nburn){
      mu_save[iter-nburn,] <- mu
      phi_save[iter-nburn,] <- phi
      sigma2_save[iter-nburn,] <- sigma2
      FF_save[iter-nburn,] <- c(FF)
      gamma_save[iter-nburn,] <- gamma
      tau_save[iter-nburn,] <- tau
    }
  }
  
  list(mu = mu_save,
       phi = phi_save,
       sigma2 = sigma2_save,
       FF = FF_save,
       gamma = gamma_save,
       tau = tau_save)
}