### Data Generation Functions ###

#' @export
make_tacfaSM_data <- function(N,
                              V,
                              M_min,
                              M_max,
                              mu,
                              sigma2,
                              rho,
                              phi_args,
                              Psi_args,
                              gamma_args,
                              tau_args){
  
  Phi <- do.call(make_Phi, phi_args)
  Psi_full <- do.call(make_Psi_null, Psi_args)
  Gamma <- do.call(make_Gamma, gamma_args)
  tau_full <- do.call(make_tau_full, tau_args)
  
  J <- ncol(Phi)
  K <- ncol(Gamma)
  L <- nrow(Phi)
  
  if (M_min == M_max){
    M <- rep(M_min, N)
  } else {
    M <- sample(M_min:M_max, N, TRUE)
  }
  
  Mu <- rep_vector(N, mu)
  R <- make_R(rho, L)
  FF <- make_FF(N, R)
  
  Y <- make_Y(Mu, FF, Phi, sigma2)
  
  Theta <- comp_Theta(tau_full, FF, Gamma)
  
  Z <- make_Z(Theta, M)
  W <- make_W(Psi_full, Z)
  
  tacfaSM_data <- list(Y = Y,
                       W = W,
                       Z = Z,
                       params = list(mu = mu,
                                     phi = phi_args$phi,
                                     FF = FF,
                                     sigma2 = sigma2,
                                     tau = tau_args$tau,
                                     gamma = gamma_args$gamma,
                                     rho = rho,
                                     Psi_full = Psi_full))
  
  tacfaSM_data
}

### Helper Functions ###

make_Psi_full <- function(Psi_full, Psi, anchor_words){
  Psi_full[-anchor_words,] <- Psi
  
  Psi_full
}

make_Phi <- function(phi, phi_factor, phi_fixed){
  L <- max(phi_factor)
  J <- length(phi_factor) + L
  
  Phi <- matrix(0, L, J)
  
  Phi[cbind(1:L, phi_fixed)] <- 1
  Phi[cbind(phi_factor, setdiff(1:J, phi_fixed))] <- phi
  
  Phi
}

make_Gamma <- function(gamma, 
                       gamma_factor,
                       gamma_fixed,
                       gamma_fixed_value,
                       topic_baseline){
  
  
  L <- max(gamma_factor)
  K <- length(gamma_factor) + L + 1
  
  Gamma <- matrix(0, L, K)
  
  Gamma[cbind(1:L, gamma_fixed)] <- gamma_fixed_value
  Gamma[cbind(gamma_factor, setdiff(1:K, c(gamma_fixed, topic_baseline)))] <- gamma
  
  Gamma
}

make_Psi_null <- function(anchor_words,
                          anchor_prob,
                          nonanchor_prob,
                          num_topic_words,
                          V = V){
  
  K <- length(anchor_words)
  Psi_null <- matrix(0, V, K)
  free_words <- setdiff(1:V, anchor_words)
  topic_word_prob <- (1-anchor_prob-(K-1)*nonanchor_prob)/(num_topic_words+1)
  
  Psi_null[anchor_words, ] <- nonanchor_prob
  Psi_null[cbind(anchor_words, 1:K)] <- anchor_prob
  
  for (k in 1:K){
    topic_words <- sample(free_words, num_topic_words)
    Psi_null[topic_words,k] <- topic_word_prob
    free_words <- setdiff(free_words, topic_words)
    Psi_null[-c(topic_words, anchor_words),k] <- (1-sum(Psi_null[,k]))/(V-K-num_topic_words)
  }
  
  Psi_null
}

make_tau_full <- function(tau, topic_baseline){
  K <- length(tau) + 1
  tau_full <- rep(0, K)
  tau_full[-topic_baseline] <- tau
  
  tau_full
}

make_R <- function(rho, L){
  R <- diag(1-rho, L) + rho
  
  R
}

make_FF <- function(N, R){
  L <- ncol(R)
  
  mu <- rep(0, L)
  
  FF <- rmvnorm(N, mu, R)
  
  FF
}

make_Y <- function(Mu, FF, Phi, sigma2){
  N <- nrow(FF)
  J <- length(sigma2)
  
  Z_mean <- numeric(J)
  Z_Sigma <- diag(sigma2)
  Z <- rmvnorm(N, Z_mean, Z_Sigma)
  
  Y <- Mu + FF%*%Phi + Z
  
  Y
}

make_W <- function(Psi, Z){
  W <- make_W_cpp(Psi, Z)
  
  W
}

make_Z <- function(Theta, M){
  Z <- make_Z_cpp(Theta, M)
  
  Z
}

comp_n_vk <- function(W, Z, V, K) {
  n_vk <- comp_n_vk_cpp(W, Z, V, K)
  
  n_vk
}

comp_n_ik <- function(Z, K){
  N <- length(Z)
  n_ik <- matrix(0, N, K)
  
  for (i in 1:N){
    n_ik[i,] <- tabulate(Z[[i]], nbins = K)
  }
  
  n_ik
}

comp_Theta <- function(tau_full, FF, Gamma){
  Theta <- comp_Theta_cpp(tau_full, FF, Gamma)
  
  Theta
}

rep_vector <- function(N, v){
  P <- length(v)
  
  vector_mat <- matrix(v, nrow = N, ncol = P, byrow = TRUE)
  
  vector_mat
}

solve_R <- function(R, chol_R) {
  backsolve(chol_R,
            forwardsolve(t(chol_R), R))
}


### Sampling Functions ###

rmvnorm <- function(n, mu, Sigma) {
  
  p <- length(mu)
  
  L <- chol(Sigma)
  
  Z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  
  samples <- Z %*% L
  
  samples <- sweep(samples, 2, mu, "+")
  
  samples
}

rtnorm_lower_single <- function(mu, sigma, lower = 0){
  repeat {
    x <- stats::rnorm(1, mu, sigma)
    
    if (x > lower) {
      return(x)
    }
  }
}

rdirichlet_single <- function(alpha){
  n <- length(alpha)
  x <- stats::rgamma(n, alpha)
  
  draw <- x/sum(x)
  
  draw
}

### Log Full Conditional Functions ###

comp_log_p_FF <- function(Y, n_ik, FF, mu, R_inv, sigma2, Phi, Gamma, tau_full, FF_Gamma){
  log_p_FF <- comp_log_p_FF_cpp(Y, n_ik, FF, mu, R_inv, sigma2, Phi, Gamma, tau_full, FF_Gamma)
  
  log_p_FF
}

comp_grad_log_p_FF <- function(Y, n_ik, FF, mu, R_inv, sigma2, Phi, Gamma, tau_full, FF_Gamma){
  grad_log_p_FF <- comp_grad_log_p_FF_cpp(Y, n_ik, FF, mu, R_inv, sigma2, Phi, 
                                          Gamma, tau_full, FF_Gamma)
  
  grad_log_p_FF
}

comp_log_p_gamma <- function(FF_Gamma, gamma, n_ik, tau, mu_gamma, sigma_gamma2){
  log_p_gamma <- comp_log_p_gamma_cpp(FF_Gamma, gamma, n_ik, tau, mu_gamma, sigma_gamma2)
  
  log_p_gamma
}

comp_grad_log_p_gamma <- function(FF,
                                  FF_Gamma,
                                  gamma,
                                  n_ik,
                                  tau_full,
                                  mu_gamma,
                                  sigma_gamma2,
                                  gamma_fixed,
                                  gamma_factor){
  
  grad_log_p_gamma <- comp_grad_log_p_gamma_cpp(FF, FF_Gamma, gamma, n_ik, tau_full, mu_gamma,
                                                sigma_gamma2, gamma_fixed, gamma_factor)
  
  grad_log_p_gamma
}

comp_log_p_tau <- function(FF_Gamma, n_ik, tau, tau_full, mu_tau, sigma_tau2){
  log_p_tau <- comp_log_p_tau_cpp(FF_Gamma, n_ik, tau, tau_full, mu_tau, sigma_tau2)
  
  log_p_tau
}

comp_grad_log_p_tau <- function(FF_Gamma, n_ik, tau, tau_full, mu_tau,
                                sigma_tau2, topic_baseline){
  grad_log_p_tau <- comp_grad_log_p_tau_cpp(FF_Gamma, n_ik, tau, tau_full, mu_tau,
                                            sigma_tau2, topic_baseline)
  
  grad_log_p_tau
}

comp_log_p_rho <- function(FF, rho, R, a_rho, b_rho) {
  S <- crossprod(FF)  # F^T F
  
  N <- nrow(FF)
  L <- ncol(FF)
  
  # Cholesky decomposition (stable + efficient)
  chol_R <- chol(R)
  
  # log |R|
  log_det_R <- 2 * sum(log(diag(chol_R)))
  
  # Compute R^{-1} S via solves (avoid explicit inverse)
  R_inv_S <- backsolve(chol_R,
                       forwardsolve(t(chol_R), S))
  
  trace_term <- sum(diag(R_inv_S))
  
  # Final log posterior (up to constant)
  log_p_rho <- -0.5 * N * log_det_R -
    0.5 * trace_term +
    (a_rho - 1) * log(rho) +
    (b_rho - 1) * log(1 - rho)
  
  log_p_rho
}

comp_grad_log_p_rho <- function(FF, rho, R, a_rho, b_rho){
  S <- crossprod(FF)
  
  N <- nrow(FF)
  L <- ncol(FF)
  
  # Build matrices
  dR <- -diag(L) + matrix(1, L, L)
  
  # Cholesky decomposition
  chol_R <- chol(R)
  
  # Compute needed pieces
  Rinv_dR <- solve_R(dR, chol_R)
  Rinv_S  <- solve_R(S, chol_R)
  
  # Trace terms
  term1 <- sum(diag(Rinv_dR))
  term2 <- sum(diag(Rinv_dR %*% Rinv_S))
  
  # Gradient
  grad_log_p_rho <- -0.5 * N * term1 +
    0.5 * term2 +
    (a_rho - 1) / rho -
    (b_rho - 1) / (1 - rho)
  
  grad_log_p_rho
}
