### Sampling method for mu ###

sample_mu <- function(Y, FF, Phi, sigma2, mu_0, sigma_mu2){
  mu <- sample_mu_cpp(Y, FF, Phi, sigma2, mu_0, sigma_mu2)
  
  mu
}

### Sampling method for sigma2

sample_sigma2 <- function(Y, FF, Phi, mu, a_sigma2, b_sigma2){
  sigma2 <- sample_sigma2_cpp(Y, FF, Phi, mu, a_sigma2, b_sigma2)
  
  sigma2
}

### Sampling method for phi ###

sample_phi <- function(Y, FF, mu, sigma2, mu_phi, sigma_phi2, phi_fixed, phi_factor){
  phi <- sample_phi_cpp(Y, FF, mu, sigma2, mu_phi, sigma_phi2, phi_fixed, phi_factor)
  
  phi
}

### Sampling method for Psi ###

sample_Psi <- function(beta,
                       n_vk,
                       anchor_words,
                       anchor_prob,
                       nonanchor_prob){
  Psi <- sample_Psi_cpp(beta, n_vk, anchor_words, anchor_prob, nonanchor_prob)
  
  Psi
}

### Sampling method for Z ###

sample_Z <- function(W, FF, tau_full, Gamma, Psi_full){
  Theta <- comp_Theta(tau_full, FF, Gamma)
  log_Psi_full <- log(Psi_full)
  
  Z <- sample_Z_cpp(W, Theta, log_Psi_full)
  
  Z
}

### Sampling method for FF ##

sample_FF <- function(
    FF, Y, n_ik, mu, R_inv, sigma2, Phi, Gamma, tau_full,
    log_step,
    iter,
    warmup,
    G_FF,
    G_inv_FF,
    FF_mean,
    C_FF,
    delta = 1e-05,
    target_accept = 0.574,
    kappa = 0.6
){
  N <- nrow(FF)
  L <- ncol(FF)
  
  I_L <- diag(L)

  step_size <- exp(log_step)
  chol_G_FF <- chol(G_FF)
  
  # ---- Current quantities ----
  FF_Gamma <- FF%*%Gamma
  
  logp_curr <- comp_log_p_FF(Y, n_ik, FF, mu, R_inv, sigma2, Phi, Gamma, tau_full, FF_Gamma)
  grad_curr <- comp_grad_log_p_FF(Y, n_ik, FF, mu, R_inv, sigma2, Phi, Gamma, tau_full, FF_Gamma)
  
  drift_curr <- grad_curr %*% G_FF
  
  # ---- Proposal ----
  noise <- matrix(stats::rnorm(N * L), N, L) %*% chol_G_FF
  
  FF_prop <- FF +
    (step_size^2 / 2) * drift_curr +
    step_size * noise
  
  # ---- Proposed quantities ----
  FF_Gamma_prop <- FF_prop%*%Gamma
  
  logp_prop <- comp_log_p_FF(Y, n_ik, FF_prop, mu, R_inv, sigma2, Phi, Gamma, tau_full, FF_Gamma_prop)
  grad_prop <- comp_grad_log_p_FF(Y, n_ik, FF_prop, mu, R_inv, sigma2, Phi, Gamma, tau_full, FF_Gamma_prop)
  
  drift_prop <- grad_prop %*% G_FF
  
  # ---- Proposal densities (Mahalanobis) ----
  mean_forward <- FF + (step_size^2 / 2) * drift_curr
  mean_backward <- FF_prop + (step_size^2 / 2) * drift_prop
  
  diff_forward <- FF_prop - mean_forward
  diff_backward <- FF - mean_backward
  
  log_q_forward <- -rowSums((diff_forward %*% G_inv_FF) * diff_forward) / (2 * step_size^2)
  log_q_backward <- -rowSums((diff_backward %*% G_inv_FF) * diff_backward) / (2 * step_size^2)
  
  # ---- MH ratio ----
  log_accept_ratio <- logp_prop - logp_curr +
    log_q_backward - log_q_forward
  
  log_u <- log(stats::runif(N))
  accept <- log_u < log_accept_ratio
  
  FF_new <- FF
  FF_new[accept, ] <- FF_prop[accept, ]
  
  acc_rate <- mean(accept)
  
  # ---- Robbins–Monro adaptation ----
  if(iter <= warmup){
    eta_t <- iter^(-kappa)
    log_step <- log_step + eta_t * (acc_rate - target_accept)
    
    # stability clamp
    log_step <- min(max(log_step, -10), 2)
    
    FF_mean_old <- FF_mean
    FF_col_mean <- colMeans(FF)
    FF_mean <- FF_mean + (FF_col_mean - FF_mean) / iter
    
    C_FF <- C_FF + tcrossprod(FF_col_mean - FF_mean_old,
                              FF_col_mean - FF_mean)
    
    if(iter > 50){
      G_FF <- C_FF / (iter - 1) + delta * I_L
      
      # Precompute inverse for next iteration
      chol_G_FF <- chol(G_FF)
      G_inv_FF <- chol2inv(chol_G_FF)
    }
  }
  
  list(
    FF = FF_new,
    log_step = log_step,
    FF_mean = FF_mean,
    C_FF = C_FF,
    G_FF = G_FF,
    G_inv_FF = G_inv_FF
  )
}

### Sampling method for gamma ###

sample_gamma <- function(
    FF, gamma, n_ik, tau_full, mu_gamma, sigma_gamma2,
    gamma_args,
    log_step,
    iter,
    warmup,
    G_gamma,
    G_inv_gamma,
    gamma_mean,
    C_gamma,
    delta = 1e-05,
    target_accept = 0.574,
    kappa = 0.6
){
  
  gamma_args$gamma <- gamma
  gamma_fixed <- gamma_args$gamma_fixed
  gamma_factor <- gamma_args$gamma_factor
  
  Gamma <- do.call(make_Gamma, gamma_args)
  
  d <- length(gamma)
  I_d <- diag(d)
  
  step_size <- exp(log_step)
  chol_G_gamma <- chol(G_gamma)
  
  # ---- Current quantities ----
  FF_Gamma <- FF%*%Gamma
  logp_curr <- comp_log_p_gamma(FF_Gamma, gamma, n_ik, tau_full, mu_gamma, sigma_gamma2)
  grad_curr <- comp_grad_log_p_gamma(FF, FF_Gamma, gamma, n_ik, tau_full, mu_gamma, sigma_gamma2, 
                                     gamma_fixed, gamma_factor)
  
  drift_curr <- grad_curr %*% G_gamma
  
  # ---- Proposal ----
  noise <- stats::rnorm(d) %*% chol_G_gamma
  
  gamma_prop <- gamma +
    (step_size^2 / 2) * drift_curr +
    step_size * noise
  
  gamma_args$gamma <- gamma_prop
  Gamma_prop <- do.call(make_Gamma, gamma_args)
  
  # ---- Proposed quantities ----
  FF_Gamma_prop <- FF%*%Gamma_prop
  logp_prop <- comp_log_p_gamma(FF_Gamma_prop, gamma_prop, n_ik, tau_full, mu_gamma, sigma_gamma2)
  grad_prop <- comp_grad_log_p_gamma(FF, FF_Gamma_prop, gamma_prop, n_ik, tau_full, mu_gamma, sigma_gamma2, 
                                     gamma_fixed, gamma_factor)
  
  drift_prop <- grad_prop %*% G_gamma
  
  # ---- Proposal densities (Mahalanobis) ----
  mean_forward <- gamma + (step_size^2 / 2) * drift_curr
  mean_backward <- gamma_prop + (step_size^2 / 2) * drift_prop
  
  diff_forward <- gamma_prop - mean_forward
  diff_backward <- gamma - mean_backward
  
  log_q_forward <- -rowSums((diff_forward %*% G_inv_gamma) * diff_forward) / (2 * step_size^2)
  log_q_backward <- -rowSums((diff_backward %*% G_inv_gamma) * diff_backward) / (2 * step_size^2)
  
  # ---- MH ratio ----
  log_accept_ratio <- logp_prop - logp_curr +
    log_q_backward - log_q_forward
  
  log_u <- log(stats::runif(1))
  accept <- log_u < log_accept_ratio
  
  if (accept){
    gamma_new <- gamma_prop
    Gamma_new <- Gamma_prop
  } else {
    gamma_new <- gamma
    Gamma_new <- Gamma
  }
  
  # ---- Robbins–Monro adaptation ----
  if(iter <= warmup){
    eta_t <- iter^(-kappa)
    log_step <- log_step + eta_t * (accept - target_accept)
    
    # stability clamp
    log_step <- min(max(log_step, -10), 2)
    
    c_gamma_new <- c(gamma_new)
    
    gamma_mean_old <- gamma_mean
    gamma_mean <- gamma_mean + (c_gamma_new - gamma_mean) / iter
    
    C_gamma <- C_gamma + tcrossprod(c_gamma_new - gamma_mean_old,
                                    c_gamma_new - gamma_mean)
    
    if(iter > 50){
      G_gamma <- C_gamma / (iter - 1) + delta * I_d
      
      # Precompute inverse for next iteration
      chol_G_gamma <- chol(G_gamma)
      G_inv_gamma <- chol2inv(chol_G_gamma)
    }
  }
  
  list(
    gamma = gamma_new,
    Gamma = Gamma_new,
    log_step = log_step,
    gamma_mean = gamma_mean,
    C_gamma = C_gamma,
    G_gamma = G_gamma,
    G_inv_gamma = G_inv_gamma
  )
}

### Sampling method for tau ###

sample_tau <- function(
    FF_Gamma, tau, n_ik, mu_tau, sigma_tau2,
    tau_args,
    log_step,
    iter,
    warmup,
    G_tau,
    G_inv_tau,
    tau_mean,
    C_tau,
    delta = 1e-05,
    target_accept = 0.574,
    kappa = 0.6
){
  
  tau_args$tau <- tau
  topic_baseline <- tau_args$topic_baseline
  
  tau_full <- do.call(make_tau_full, tau_args)
  
  d <- length(tau)
  I_d <- diag(d)
  
  step_size <- exp(log_step)
  chol_G_tau <- chol(G_tau)
  
  # ---- Current quantities ----
  logp_curr <- comp_log_p_tau(FF_Gamma, n_ik, tau, tau_full, mu_tau, sigma_tau2)
  grad_curr <- comp_grad_log_p_tau(FF_Gamma, n_ik, tau, tau_full, mu_tau, sigma_tau2, topic_baseline)
  
  drift_curr <- grad_curr %*% G_tau
  
  # ---- Proposal ----
  noise <- stats::rnorm(d) %*% chol_G_tau
  
  tau_prop <- tau +
    (step_size^2 / 2) * drift_curr +
    step_size * noise
  
  tau_args$tau <- tau_prop
  tau_full_prop <- do.call(make_tau_full, tau_args)
  
  # ---- Proposed quantities ----
  logp_prop <-  comp_log_p_tau(FF_Gamma, n_ik, tau_prop, tau_full_prop, mu_tau, sigma_tau2)
  grad_prop <- comp_grad_log_p_tau(FF_Gamma, n_ik, tau_prop, tau_full_prop, mu_tau, sigma_tau2, topic_baseline)
  
  drift_prop <- grad_prop %*% G_tau
  
  # ---- Proposal densities (Mahalanobis) ----
  mean_forward <- tau + (step_size^2 / 2) * drift_curr
  mean_backward <- tau_prop + (step_size^2 / 2) * drift_prop
  
  diff_forward <- tau_prop - mean_forward
  diff_backward <- tau - mean_backward
  
  log_q_forward <- -rowSums((diff_forward %*% G_inv_tau) * diff_forward) / (2 * step_size^2)
  log_q_backward <- -rowSums((diff_backward %*% G_inv_tau) * diff_backward) / (2 * step_size^2)
  
  # ---- MH ratio ----
  log_accept_ratio <- logp_prop - logp_curr +
    log_q_backward - log_q_forward
  
  log_u <- log(stats::runif(1))
  accept <- log_u < log_accept_ratio
  
  if (accept){
    tau_new <- tau_prop
    tau_full_new <- tau_full_prop
  } else {
    tau_new <- tau
    tau_full_new <- tau_full
  }
  
  # ---- Robbins–Monro adaptation ----
  if(iter <= warmup){
    eta_t <- iter^(-kappa)
    log_step <- log_step + eta_t * (accept - target_accept)
    
    # stability clamp
    log_step <- min(max(log_step, -10), 2)
    
    c_tau_new <- c(tau_new)
    
    tau_mean_old <- tau_mean
    tau_mean <- tau_mean + (c_tau_new - tau_mean) / iter
    
    C_tau <- C_tau + tcrossprod(c_tau_new - tau_mean_old,
                                c_tau_new - tau_mean)
    
    if(iter > 50){
      G_tau <- C_tau / (iter - 1) + delta * I_d
      
      # Precompute inverse for next iteration
      chol_G_tau <- chol(G_tau)
      G_inv_tau <- chol2inv(chol_G_tau)
    }
  }
  
  list(
    tau = tau_new,
    tau_full = tau_full_new,
    log_step = log_step,
    tau_mean = tau_mean,
    C_tau = C_tau,
    G_tau = G_tau,
    G_inv_tau = G_inv_tau
  )
}

### Sampling method for rho ###

sample_rho <- function(
    rho,
    R,
    FF,
    a_rho,
    b_rho,
    log_step,
    iter,
    warmup,
    target_accept = 0.574,
    kappa = 0.6
){
  
  L <- ncol(FF)
  
  # ---- Transform ----
  theta <- stats::qlogis(rho)
  step_size <- exp(log_step)
  
  logp_curr <- comp_log_p_rho(FF, rho, R, a_rho, b_rho) +
    log(rho * (1 - rho))  # Jacobian
  
  grad_rho_curr <- comp_grad_log_p_rho(FF, rho, R, a_rho, b_rho)
  
  grad_theta_curr <- grad_rho_curr * rho * (1 - rho) +
    (1 - 2 * rho)
  
  drift_curr <- grad_theta_curr
  
  # ---- Proposal ----
  noise <- rnorm(1)
  
  theta_prop <- theta +
    (step_size^2 / 2) * drift_curr +
    step_size * noise
  
  rho_prop <- stats::plogis(theta_prop)
  
  # ---- Proposed quantities ----
  R_prop <- make_R(rho_prop, L)
  
  logp_prop <- comp_log_p_rho(FF, rho_prop, R_prop, a_rho, b_rho) +
    log(rho_prop * (1 - rho_prop))
  
  grad_rho_prop <- comp_grad_log_p_rho(FF, rho_prop, R_prop, a_rho, b_rho)
  
  grad_theta_prop <- grad_rho_prop * rho_prop * (1 - rho_prop) +
    (1 - 2 * rho_prop)
  
  drift_prop <- grad_theta_prop
  
  # ---- Proposal densities ----
  mean_forward <- theta + (step_size^2 / 2) * drift_curr
  mean_backward <- theta_prop + (step_size^2 / 2) * drift_prop
  
  diff_forward <- theta_prop - mean_forward
  diff_backward <- theta - mean_backward
  
  log_q_forward <- - (diff_forward^2) / (2 * step_size^2)
  log_q_backward <- - (diff_backward^2) / (2 * step_size^2)
  
  # ---- MH ratio ----
  log_accept_ratio <- logp_prop - logp_curr +
    log_q_backward - log_q_forward
  
  accept <- log(runif(1)) < log_accept_ratio
  
  if (accept){
    rho_new <- rho_prop
    theta <- theta_prop
    R_new <- R_prop
  } else {
    rho_new <- rho
    R_new <- R
  }
  
  # ---- Adaptation ----
  if(iter <= warmup){
    eta_t <- iter^(-kappa)
    log_step <- log_step + eta_t * (accept - target_accept)
    log_step <- min(max(log_step, -10), 2)
  }
  
  list(
    rho = rho_new,
    R = R_new,
    log_step = log_step
  )
}
