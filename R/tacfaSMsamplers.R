sample_z <- function(
    W, z,
    F_mat, tau, gamma,
    psi_star,
    lambda_T,
    B, c_val,
    n_vk, n_ik
) {
  
  N <- length(W)
  K <- length(tau)
  
  log_psi_star <- log(psi_star)
  
  for (i in seq_len(N)) {
    
    words <- W[[i]]
    z_i   <- z[[i]]
    M_i   <- length(words)
    
    f_i <- F_mat[i, ]
    
    ## ----------------------------------------
    ## Compute log theta_i (once per document)
    ## ----------------------------------------
    
    tau_star <- tau
    tau_star[K] <- 0
    
    gamma_star <- gamma
    gamma_star[K] <- 0
    gamma_star[B] <- c_val
    
    log_theta <- tau_star + gamma_star * f_i[lambda_T]
    log_theta <- log_theta - log_sum_exp(log_theta)
    
    log_p <- numeric(K)
    
    ## ----------------------------------------
    ## Loop over tokens
    ## ----------------------------------------
    
    for (m in seq_len(M_i)) {
      
      v <- words[m]
      k_old <- z_i[m]
      
      ## remove old counts
      n_vk[v, k_old] <- n_vk[v, k_old] - 1
      n_ik[i, k_old] <- n_ik[i, k_old] - 1
      
      
      ## ------------------------------------
      ## Compute log probabilities
      ## ------------------------------------
      
      log_p[] <- log_theta + log_psi_star[v, ]
      
      # numerically stable normalization
      log_p <- log_p - max(log_p)
      p <- exp(log_p)
      p <- p / sum(p)
      
      
      ## ------------------------------------
      ## Sample new topic
      ## ------------------------------------
      
      k_new <- sample_cat(p)
      
      z_i[m] <- k_new
      
      
      ## ------------------------------------
      ## update counts
      ## ------------------------------------
      
      n_vk[v, k_new] <- n_vk[v, k_new] + 1
      n_ik[i, k_new] <- n_ik[i, k_new] + 1
    }
    
    z[[i]] <- z_i
  }
  
  
  return(list(
    z = z,
    n_vk = n_vk,
    n_ik = n_ik
  ))
}

sample_psi <- function(
    n_vk,
    V_set,
    hyper
) {
  
  beta    <- hyper$beta
  eta     <- hyper$eta
  epsilon <- hyper$epsilon
  
  V <- nrow(n_vk)
  K <- ncol(n_vk)
  
  remaining_mass <- 1 - eta - (K - 1) * epsilon
  
  
  ### ------------------------------------------------------------
  ### Setup
  ### ------------------------------------------------------------
  
  psi      <- matrix(0, V, K)
  psi_star <- matrix(0, V, K)
  
  non_V <- setdiff(seq_len(V), V_set)
  
  
  ### ------------------------------------------------------------
  ### Loop over topics
  ### ------------------------------------------------------------
  
  for (k in seq_len(K)) {
    
    v_k <- V_set[k]
    
    ## ----------------------------------------
    ## 1. Fixed anchor structure
    ## ----------------------------------------
    
    psi_star[v_k, k] <- eta
    
    other_V <- setdiff(V_set, v_k)
    psi_star[other_V, k] <- epsilon
    
    
    ## ----------------------------------------
    ## 2. Dirichlet update for non-V
    ## ----------------------------------------
    
    alpha_post <- beta[non_V] + n_vk[non_V, k]
    
    dir_draw <- rdirichlet(alpha_post)
    
    psi[non_V, k] <- remaining_mass * dir_draw
    
    psi_star[non_V, k] <- psi[non_V, k]
  }
  
  
  ### ------------------------------------------------------------
  ### Return
  ### ------------------------------------------------------------
  
  return(list(
    psi = psi,
    psi_star = psi_star
  ))
}

sample_mu <- function(
    Y, F_mat,
    phi,
    sigma2,
    lambda_S,
    A,
    hyper
) {
  
  mu_0      <- hyper$mu_0
  sigma_mu2 <- hyper$sigma_mu2
  
  N <- nrow(Y)
  J <- ncol(Y)
  
  
  ### ------------------------------------------------------------
  ### 1. Compute phi*
  ### ------------------------------------------------------------
  
  phi_star <- phi
  phi_star[A] <- 1
  
  
  ### ------------------------------------------------------------
  ### 2. Compute residual sums
  ### ------------------------------------------------------------
  
  # residual sum for each j
  sum_resid <- numeric(J)
  
  for (j in seq_len(J)) {
    
    l_j <- lambda_S[j]
    
    # vector over i
    resid_ij <- Y[, j] - phi_star[j] * F_mat[, l_j]
    
    sum_resid[j] <- sum(resid_ij)
  }
  
  
  ### ------------------------------------------------------------
  ### 3. Posterior parameters
  ### ------------------------------------------------------------
  
  precision <- N / sigma2 + 1 / sigma_mu2
  var_post  <- 1 / precision
  
  mean_post <- (
    (sum_resid / sigma2) +
      (mu_0 / sigma_mu2)
  ) / precision
  
  
  ### ------------------------------------------------------------
  ### 4. Sample
  ### ------------------------------------------------------------
  
  mu <- rnorm(J, mean = mean_post, sd = sqrt(var_post))
  
  
  return(mu)
}

sample_sigma2 <- function(
    Y, F_mat,
    mu, phi,
    lambda_S,
    A,
    hyper
) {
  
  a_sigma2 <- hyper$a_sigma2
  b_sigma2 <- hyper$b_sigma2
  
  N <- nrow(Y)
  J <- ncol(Y)
  
  
  ### ------------------------------------------------------------
  ### 1. Compute phi*
  ### ------------------------------------------------------------
  
  phi_star <- phi
  phi_star[A] <- 1
  
  
  ### ------------------------------------------------------------
  ### 2. Compute SSE for each j
  ### ------------------------------------------------------------
  
  sse <- numeric(J)
  
  for (j in seq_len(J)) {
    
    l_j <- lambda_S[j]
    
    resid_ij <- Y[, j] - mu[j] - phi_star[j] * F_mat[, l_j]
    
    sse[j] <- sum(resid_ij^2)
  }
  
  
  ### ------------------------------------------------------------
  ### 3. Posterior parameters
  ### ------------------------------------------------------------
  
  shape_post <- a_sigma2 + N / 2
  rate_post  <- b_sigma2 + 0.5 * sse
  
  
  ### ------------------------------------------------------------
  ### 4. Sample (InvGamma via Gamma)
  ### ------------------------------------------------------------
  
  sigma2 <- 1 / rgamma(J, shape = shape_post, rate = rate_post)
  
  
  return(sigma2)
}

sample_phi <- function(
    Y, F_mat,
    mu, sigma2,
    lambda_S,
    A,
    hyper
) {
  
  mu_phi     <- hyper$mu_phi
  sigma_phi2 <- hyper$sigma_phi2
  
  N <- nrow(Y)
  J <- ncol(Y)
  
  
  ### ------------------------------------------------------------
  ### Initialize
  ### ------------------------------------------------------------
  
  phi <- numeric(J)
  
  # fixed indices
  phi[A] <- 1
  
  # free indices
  free_idx <- setdiff(seq_len(J), A)
  
  
  ### ------------------------------------------------------------
  ### Loop over free j
  ### ------------------------------------------------------------
  
  for (j in free_idx) {
    
    l_j <- lambda_S[j]
    
    f_j <- F_mat[, l_j]  # vector length N
    
    ## ----------------------------------------
    ## Compute sums
    ## ----------------------------------------
    
    y_centered <- Y[, j] - mu[j]
    
    sum_fy <- sum(f_j * y_centered)
    sum_f2 <- sum(f_j^2)
    
    
    ## ----------------------------------------
    ## Posterior parameters
    ## ----------------------------------------
    
    precision <- (sum_f2 / sigma2[j]) + (1 / sigma_phi2)
    var_post  <- 1 / precision
    
    mean_post <- (
      (sum_fy / sigma2[j]) +
        (mu_phi / sigma_phi2)
    ) / precision
    
    
    ## ----------------------------------------
    ## Sample truncated normal
    ## ----------------------------------------
    
    phi[j] <- rtruncnorm_pos(
      n = 1,
      mean = mean_post,
      sd = sqrt(var_post)
    )
  }
  
  phi[A] <- 1
  
  return(phi)
}

sample_F <- function(
    F_mat,
    adapt_f,
    t_global,
    T_adapt,
    a,
    Y,
    mu, phi, sigma2,
    lambda_S,
    z,
    tau, gamma,
    lambda_T,
    B, c_val,
    R_inv,
    A,
    target_accept = 0.3
) {
  
  N <- nrow(F_mat)
  L <- ncol(F_mat)
  d <- L
  
  adapt_phase <- (t_global <= T_adapt)
  
  accept_vec <- numeric(N)
  
  tau_star <- tau
  tau_star[length(tau)] <- 0
  
  gamma_star <- gamma
  gamma_star[length(gamma)] <- 0
  gamma_star[B] <- c_val
  
  ### ------------------------------------------------------------
  ### Loop over i
  ### ------------------------------------------------------------
  
  for (i in seq_len(N)) {
    
    f_i <- F_mat[i, ]
    
    state_i <- adapt_f[[i]]
    
    log_s <- state_i$log_s
    s     <- exp(log_s)
    
    mean_i <- state_i$mean
    cov_i  <- state_i$cov
    t_i    <- state_i$t
    
    
    ## ----------------------------------------------------------
    ## Build proposal covariance
    ## ----------------------------------------------------------
    
    if (t_i > 1) {
      cov_est <- cov_i / (t_i - 1)
    } else {
      cov_est <- diag(L)
    }
    
    Sigma <- (2.38^2 / d) * build_proposal_cov(s, cov_est)
    
    
    ## ----------------------------------------------------------
    ## Proposal
    ## ----------------------------------------------------------
    
    L_Sigma <- chol(Sigma)
    eps <- rnorm(d)
    f_prop <- f_i + as.numeric(L_Sigma %*% eps)
    
    
    ## ----------------------------------------------------------
    ## Log posterior
    ## ----------------------------------------------------------
    
    log_curr <- log_p_f_i(
      f_i = f_i,
      i = i,
      Y = Y,
      mu = mu,
      phi = phi,
      sigma2 = sigma2,
      lambda_S = lambda_S,
      z_i = z[[i]],
      tau_star = tau_star,
      gamma_star = gamma_star,
      lambda_T = lambda_T,
      B = B,
      c_val = c_val,
      R_inv = R_inv,
      A = A
    )
    
    log_prop <- log_p_f_i(
      f_i = f_prop,
      i = i,
      Y = Y,
      mu = mu,
      phi = phi,
      sigma2 = sigma2,
      lambda_S = lambda_S,
      z_i = z[[i]],
      tau_star = tau_star,
      gamma_star = gamma_star,
      lambda_T = lambda_T,
      B = B,
      c_val = c_val,
      R_inv = R_inv,
      A = A
    )
    
    
    ## ----------------------------------------------------------
    ## Accept / Reject
    ## ----------------------------------------------------------
    
    log_alpha <- log_prop - log_curr
    accept <- (log(runif(1)) < log_alpha)
    
    if (accept) {
      F_mat[i, ] <- f_prop
      f_i_new <- f_prop
      accept_vec[i] <- 1
    } else {
      f_i_new <- f_i
      accept_vec[i] <- 0
    }
    
    
    ## ----------------------------------------------------------
    ## Adaptation (only during burn-in)
    ## ----------------------------------------------------------
    
    if (adapt_phase) {
      
      ## Robbins–Monro step size
      kappa_t <- t_global^(-a)
      
      ## update scale
      log_s_new <- rm_update(
        log_s,
        alpha = accept_vec[i],
        target = target_accept,
        kappa_t = kappa_t
      )
      
      ## update empirical covariance
      ecov <- empirical_cov_update(
        mean = mean_i,
        cov  = cov_i,
        x    = f_i_new,
        t    = t_i + 1
      )
      
      adapt_f[[i]] <- list(
        log_s = log_s_new,
        mean  = ecov$mean,
        cov   = ecov$cov,
        t     = t_i + 1
      )
      
    }
  }
  
  
  ### ------------------------------------------------------------
  ### Return
  ### ------------------------------------------------------------
  
  return(list(
    F_mat = F_mat,
    adapt_f = adapt_f,
    accept_rate = mean(accept_vec)
  ))
}

sample_tau <- function(
    tau,
    adapt_tau,
    t_global,
    T_adapt,
    a,
    F_mat,
    z,
    gamma,
    lambda_T,
    B, c_val,
    hyper,
    target_accept = 0.234
) {
  
  K <- length(tau)
  
  adapt_phase <- (t_global <= T_adapt)
  
  
  ### ------------------------------------------------------------
  ### Extract adaptive state
  ### ------------------------------------------------------------
  
  log_s <- adapt_tau$log_s
  s     <- exp(log_s)
  
  mean_tau <- adapt_tau$mean
  cov_tau  <- adapt_tau$cov
  t_tau    <- adapt_tau$t
  
  
  ### ------------------------------------------------------------
  ### Proposal covariance
  ### ------------------------------------------------------------
  
  d <- K - 1
  
  if (t_tau > 1) {
    cov_est <- cov_tau / (t_tau - 1)
  } else {
    cov_est <- diag(d)
  }
  
  Sigma <- (2.38^2 / d) * build_proposal_cov(s, cov_est)
  
  
  ### ------------------------------------------------------------
  ### Proposal (only first K-1)
  ### ------------------------------------------------------------
  
  tau_prop <- tau
  tau_prop[1:(K-1)] <- as.numeric(
    MASS::mvrnorm(1, mu = tau[1:(K-1)], Sigma = Sigma)
  )
  
  tau_prop[K] <- 0  # enforce baseline
  
  
  ### ------------------------------------------------------------
  ### Log posterior
  ### ------------------------------------------------------------
  
  log_curr <- log_p_tau(
    tau, F_mat, z,
    gamma, lambda_T,
    B, c_val,
    hyper
  )
  
  log_prop <- log_p_tau(
    tau_prop, F_mat, z,
    gamma, lambda_T,
    B, c_val,
    hyper
  )
  
  
  ### ------------------------------------------------------------
  ### Accept / Reject
  ### ------------------------------------------------------------
  
  log_alpha <- log_prop - log_curr
  accept <- (log(runif(1)) < log_alpha)
  
  if (accept) {
    tau_new <- tau_prop
    alpha_val <- 1
  } else {
    tau_new <- tau
    alpha_val <- 0
  }
  
  
  ### ------------------------------------------------------------
  ### Adaptation
  ### ------------------------------------------------------------
  
  if (adapt_phase) {
    
    kappa_t <- t_global^(-a)
    
    log_s_new <- rm_update(
      log_s,
      alpha = alpha_val,
      target = target_accept,
      kappa_t = kappa_t
    )
    
    ecov <- empirical_cov_update(
      mean = mean_tau,
      cov  = cov_tau,
      x    = tau_new[1:(K-1)],
      t    = t_tau + 1
    )
    
    adapt_tau <- list(
      log_s = log_s_new,
      mean  = ecov$mean,
      cov   = ecov$cov,
      t     = t_tau + 1
    )
  }
  
  
  return(list(
    tau = tau_new,
    adapt_tau = adapt_tau,
    accept = alpha_val
  ))
}

sample_gamma <- function(
    gamma,
    adapt_gamma,
    t_global,
    T_adapt,
    a,
    F_mat,
    z,
    tau,
    lambda_T,
    B, c_val,
    hyper,
    target_accept = 0.234
) {
  
  K <- length(gamma)
  
  free_idx <- setdiff(1:(K-1), B)
  d <- length(free_idx)
  
  adapt_phase <- (t_global <= T_adapt)
  
  
  ### ------------------------------------------------------------
  ### Extract adaptive state
  ### ------------------------------------------------------------
  
  log_s <- adapt_gamma$log_s
  s     <- exp(log_s)
  
  mean_g <- adapt_gamma$mean
  cov_g  <- adapt_gamma$cov
  t_g    <- adapt_gamma$t
  
  
  ### ------------------------------------------------------------
  ### Proposal covariance
  ### ------------------------------------------------------------
  
  if (t_g > 1) {
    cov_est <- cov_g / (t_g - 1)
  } else {
    cov_est <- diag(d)
  }
  
  if (d > 0) {
    Sigma <- (2.38^2 / d) * build_proposal_cov(s, cov_est)
  }  
  
  ### ------------------------------------------------------------
  ### Proposal (only free indices)
  ### ------------------------------------------------------------
  
  gamma_prop <- gamma
  
  gamma_prop[free_idx] <- as.numeric(
    MASS::mvrnorm(1, mu = gamma[free_idx], Sigma = Sigma)
  )
  
  # enforce constraints
  gamma_prop[B] <- c_val
  gamma_prop[K] <- 0
  
  
  ### ------------------------------------------------------------
  ### Log posterior
  ### ------------------------------------------------------------
  
  log_curr <- log_p_gamma(
    gamma,
    F_mat,
    z,
    tau,
    lambda_T,
    B, c_val,
    hyper
  )
  
  log_prop <- log_p_gamma(
    gamma_prop,
    F_mat,
    z,
    tau,
    lambda_T,
    B, c_val,
    hyper
  )
  
  
  ### ------------------------------------------------------------
  ### Accept / Reject
  ### ------------------------------------------------------------
  
  log_alpha <- log_prop - log_curr
  accept <- (log(runif(1)) < log_alpha)
  
  if (accept) {
    gamma_new <- gamma_prop
    alpha_val <- 1
  } else {
    gamma_new <- gamma
    alpha_val <- 0
  }
  
  
  ### ------------------------------------------------------------
  ### Adaptation
  ### ------------------------------------------------------------
  
  if (adapt_phase && d > 0) {
    
    kappa_t <- t_global^(-a)
    
    log_s_new <- rm_update(
      log_s,
      alpha = alpha_val,
      target = target_accept,
      kappa_t = kappa_t
    )
    
    ecov <- empirical_cov_update(
      mean = mean_g,
      cov  = cov_g,
      x    = gamma_new[free_idx],
      t    = t_g + 1
    )
    
    adapt_gamma <- list(
      log_s = log_s_new,
      mean  = ecov$mean,
      cov   = ecov$cov,
      t     = t_g + 1
    )
  }
  
  
  return(list(
    gamma = gamma_new,
    adapt_gamma = adapt_gamma,
    accept = alpha_val
  ))
}

sample_rho <- function(
    rho,
    adapt_rho,
    t_global,
    T_adapt,
    a,
    F_mat,
    hyper,
    L,
    target_accept = 0.44
) {
  
  adapt_phase <- (t_global <= T_adapt)
  
  
  ### ------------------------------------------------------------
  ### Transform
  ### ------------------------------------------------------------
  
  xi <- log(rho / (1 - rho))
  
  
  ### ------------------------------------------------------------
  ### Proposal
  ### ------------------------------------------------------------
  
  log_s <- adapt_rho$log_s
  s     <- exp(log_s)
  
  xi_prop <- rnorm(1, mean = xi, sd = s)
  
  rho_prop <- 1 / (1 + exp(-xi_prop))
  
  
  ### ------------------------------------------------------------
  ### Log posterior (with Jacobian)
  ### ------------------------------------------------------------
  
  log_curr <- log_p_rho(rho, F_mat, hyper, L) +
    log(rho) + log(1 - rho)
  
  log_prop <- log_p_rho(rho_prop, F_mat, hyper, L) +
    log(rho_prop) + log(1 - rho_prop)
  
  
  ### ------------------------------------------------------------
  ### Accept / Reject
  ### ------------------------------------------------------------
  
  log_alpha <- log_prop - log_curr
  accept <- (log(runif(1)) < log_alpha)
  
  if (accept) {
    rho_new <- rho_prop
    xi_new  <- xi_prop
    alpha_val <- 1
  } else {
    rho_new <- rho
    xi_new  <- xi
    alpha_val <- 0
  }
  
  
  ### ------------------------------------------------------------
  ### Adaptation (scalar)
  ### ------------------------------------------------------------
  
  if (adapt_phase) {
    
    kappa_t <- t_global^(-a)
    
    log_s_new <- rm_update(
      log_s,
      alpha = alpha_val,
      target = target_accept,
      kappa_t = kappa_t
    )
    
    adapt_rho <- list(
      log_s = log_s_new
    )
  }
  
  
  ### ------------------------------------------------------------
  ### Recompute R
  ### ------------------------------------------------------------
  
  Rstuff <- compute_R(rho_new, L)
  
  
  return(list(
    rho = rho_new,
    adapt_rho = adapt_rho,
    R = Rstuff$R,
    R_inv = Rstuff$R_inv,
    log_det_R = Rstuff$log_det_R,
    accept = alpha_val
  ))
}