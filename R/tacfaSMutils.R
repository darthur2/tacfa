sample_cat <- function(p){
  sample_catC(p)
}

flatten_param_chains <- function(arr, name) {
  
  dims <- dim(arr)
  
  # Expected: (iter, chain, ...)
  n_iter  <- dims[1]
  n_chain <- dims[2]
  param_dims <- dims[-c(1,2)]
  
  n_param <- prod(param_dims)
  
  mat <- array(arr, dim = c(n_iter, n_chain, n_param))
  
  # build parameter names
  index_grid <- expand.grid(lapply(param_dims, seq_len))
  
  param_names <- apply(index_grid, 1, function(idx) {
    paste0(name, "[", paste(idx, collapse = ","), "]")
  })
  
  dimnames(mat) <- list(NULL, NULL, param_names)
  
  return(mat)
}

build_draws_chains <- function(mu, sigma2, phi, F, tau, gamma, rho) {
  
  rho_arr <- array(rho, dim = c(dim(mu)[1], dim(mu)[2], 1))
  
  arrays <- list(
    flatten_param_chains(mu, "mu"),
    flatten_param_chains(sigma2, "sigma2"),
    flatten_param_chains(phi, "phi"),
    flatten_param_chains(F, "F"),
    flatten_param_chains(tau, "tau"),
    flatten_param_chains(gamma, "gamma"),
    flatten_param_chains(rho_arr, "rho")
  )
  
  # combine along variable dimension
  combined <- do.call(abind::abind, c(arrays, along = 3))
  
  draws <- posterior::as_draws_array(combined)
  
  return(draws)
}

rtruncnorm_pos <- function(n, mean, sd) {
  x <- rnorm(n, mean, sd)
  
  # rejection sampling
  while (any(x <= 0)) {
    idx <- which(x <= 0)
    x[idx] <- rnorm(length(idx), mean, sd)
  }
  
  return(x)
}

init_cfa_params <- function(J, A, hyper) {
  
  # unpack hyperparameters
  mu_0       <- hyper$mu_0
  sigma_mu2  <- hyper$sigma_mu2
  
  mu_phi     <- hyper$mu_phi
  sigma_phi2 <- hyper$sigma_phi2
  
  a_sigma2   <- hyper$a_sigma2
  b_sigma2   <- hyper$b_sigma2
  
  
  ### ------------------------------------------------------------
  ### 1. Sample mu_j
  ### ------------------------------------------------------------
  
  mu <- rnorm(J, mean = mu_0, sd = sqrt(sigma_mu2))
  
  
  ### ------------------------------------------------------------
  ### 2. Sample sigma_j^2
  ### ------------------------------------------------------------
  
  sigma2 <- 1 / rgamma(J, shape = a_sigma2, rate = b_sigma2)
  
  
  ### ------------------------------------------------------------
  ### 3. Sample phi_j
  ### ------------------------------------------------------------
  
  phi <- numeric(J)
  
  # indices not in A
  not_A <- setdiff(seq_len(J), A)
  
  # truncated normal for free loadings
  phi[not_A] <- rtruncnorm_pos(
    n = length(not_A),
    mean = mu_phi,
    sd = sqrt(sigma_phi2)
  )
  
  # fixed loadings
  phi[A] <- 1
  
  
  ### ------------------------------------------------------------
  ### Return
  ### ------------------------------------------------------------
  
  return(list(
    mu = mu,
    sigma2 = sigma2,
    phi = phi
  ))
}

init_factors <- function(N, L) {
  
  F_mat <- matrix(
    rnorm(N * L, mean = 0, sd = 1),
    nrow = N,
    ncol = L
  )
  
  return(F_mat)
}

init_topic_params <- function(K, B, hyper) {
  
  # unpack hyperparameters
  mu_tau      <- hyper$mu_tau
  sigma_tau2  <- hyper$sigma_tau2
  
  mu_gamma     <- hyper$mu_gamma
  sigma_gamma2 <- hyper$sigma_gamma2
  
  c_val <- hyper$c
  
  
  ### ------------------------------------------------------------
  ### 1. Initialize tau
  ### ------------------------------------------------------------
  
  tau <- numeric(K)
  
  # sample for k = 1,...,K-1
  tau[1:(K-1)] <- rnorm(
    K - 1,
    mean = mu_tau,
    sd   = sqrt(sigma_tau2)
  )
  
  # baseline
  tau[K] <- 0
  
  
  ### ------------------------------------------------------------
  ### 2. Initialize gamma
  ### ------------------------------------------------------------
  
  gamma <- numeric(K)
  
  # indices excluding baseline
  non_baseline <- 1:(K-1)
  
  # indices in B but not baseline
  B_nb <- intersect(B, non_baseline)
  
  # indices free (not in B and not baseline)
  free_idx <- setdiff(non_baseline, B_nb)
  
  # fixed values for B
  gamma[B_nb] <- c_val
  
  # sampled values for free indices
  gamma[free_idx] <- rnorm(
    length(free_idx),
    mean = mu_gamma,
    sd   = sqrt(sigma_gamma2)
  )
  
  # baseline
  gamma[K] <- 0
  
  
  ### ------------------------------------------------------------
  ### Return
  ### ------------------------------------------------------------
  
  return(list(
    tau = tau,
    gamma = gamma
  ))
}

rdirichlet <- function(alpha) {
  x <- rgamma(length(alpha), shape = alpha, rate = 1)
  x / sum(x)
}

init_psi <- function(K, V, V_set, hyper) {
  
  beta   <- hyper$beta
  eta    <- hyper$eta
  epsilon <- hyper$epsilon
  
  
  ### ------------------------------------------------------------
  ### 0. Checks
  ### ------------------------------------------------------------
  
  if (length(V_set) != K)
    stop("V_set must have length K")
  
  if (length(beta) != V)
    stop("beta must have length V")
  
  remaining_mass <- 1 - eta - (K - 1) * epsilon
  
  if (remaining_mass <= 0)
    stop("Invalid eta/epsilon: remaining mass must be positive")
  
  
  ### ------------------------------------------------------------
  ### 1. Setup
  ### ------------------------------------------------------------
  
  psi      <- matrix(0, V, K)
  psi_star <- matrix(0, V, K)
  
  non_V <- setdiff(seq_len(V), V_set)
  
  
  ### ------------------------------------------------------------
  ### 2. Loop over topics
  ### ------------------------------------------------------------
  
  for (k in seq_len(K)) {
    
    v_k <- V_set[k]
    
    ## ----------------------------------------
    ## 2.1 Anchor and special vocab
    ## ----------------------------------------
    
    psi_star[v_k, k] <- eta
    
    other_V <- setdiff(V_set, v_k)
    psi_star[other_V, k] <- epsilon
    
    
    ## ----------------------------------------
    ## 2.2 Dirichlet over non-V
    ## ----------------------------------------
    
    alpha <- beta[non_V]
    
    dir_draw <- rdirichlet(alpha)
    
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

init_assignments <- function(W, N, K, V) {
  
  ### ------------------------------------------------------------
  ### 1. Initialize z
  ### ------------------------------------------------------------
  
  z <- vector("list", N)
  
  # counts
  n_vk <- matrix(0, nrow = V, ncol = K)
  n_ik <- matrix(0, nrow = N, ncol = K)
  
  
  ### ------------------------------------------------------------
  ### 2. Loop over documents
  ### ------------------------------------------------------------
  
  for (i in seq_len(N)) {
    
    words <- W[[i]]
    M_i   <- length(words)
    
    ## ----------------------------------------
    ## Sample z_im ~ Uniform(1,...,K)
    ## ----------------------------------------
    
    z_i <- sample.int(K, M_i, replace = TRUE)
    
    z[[i]] <- z_i
    
    
    ## ----------------------------------------
    ## Update n_ik (document-topic counts)
    ## ----------------------------------------
    
    n_ik[i, ] <- tabulate(z_i, nbins = K)
    
    
    ## ----------------------------------------
    ## Update n_vk (word-topic counts)
    ## ----------------------------------------
    
    # loop over tokens (cannot avoid because words differ)
    for (m in seq_len(M_i)) {
      v <- words[m]
      k <- z_i[m]
      
      n_vk[v, k] <- n_vk[v, k] + 1
    }
  }
  
  
  ### ------------------------------------------------------------
  ### Return
  ### ------------------------------------------------------------
  
  return(list(
    z = z,
    n_vk = n_vk,
    n_ik = n_ik
  ))
}

compute_R <- function(rho, L) {
  
  ### ------------------------------------------------------------
  ### 0. Basic checks
  ### ------------------------------------------------------------
  
  if (rho <= 0 || rho >= 1) {
    warning("rho is outside (0,1); matrix may not be valid")
  }
  
  
  ### ------------------------------------------------------------
  ### 1. Construct R
  ### ------------------------------------------------------------
  
  R <- matrix(rho, nrow = L, ncol = L)
  diag(R) <- 1
  
  
  ### ------------------------------------------------------------
  ### 2. Log determinant (closed form)
  ### ------------------------------------------------------------
  
  log_det_R <- (L - 1) * log(1 - rho) +
    log(1 + (L - 1) * rho)
  
  
  ### ------------------------------------------------------------
  ### 3. Inverse (closed form)
  ### ------------------------------------------------------------
  
  one_vec <- rep(1, L)
  
  denom <- 1 + (L - 1) * rho
  
  R_inv <- (1 / (1 - rho)) * (
    diag(L) -
      (rho / denom) * (one_vec %*% t(one_vec))
  )
  
  
  ### ------------------------------------------------------------
  ### Return
  ### ------------------------------------------------------------
  
  return(list(
    R = R,
    R_inv = R_inv,
    log_det_R = log_det_R
  ))
}

initialize_parameters <- function(
    Y, W,
    lambda_S, lambda_T,
    A, B, V_set,
    S_list, T_list,
    K, L, V,
    hyper
) {
  
  ### ------------------------------------------------------------
  ### 1. Dimensions
  ### ------------------------------------------------------------
  
  N <- nrow(Y)
  J <- ncol(Y)
  
  
  ### ------------------------------------------------------------
  ### 2. CFA parameters
  ### ------------------------------------------------------------
  
  cfa <- init_cfa_params(
    J = J,
    A = A,
    hyper = hyper
  )
  
  mu     <- cfa$mu
  sigma2 <- cfa$sigma2
  phi    <- cfa$phi
  
  
  ### ------------------------------------------------------------
  ### 3. Latent factors
  ### ------------------------------------------------------------
  
  F_mat <- init_factors(N = N, L = L)
  
  
  ### ------------------------------------------------------------
  ### 4. Topic parameters
  ### ------------------------------------------------------------
  
  topic <- init_topic_params(
    K = K,
    B = B,
    hyper = hyper
  )
  
  tau   <- topic$tau
  gamma <- topic$gamma
  
  
  ### ------------------------------------------------------------
  ### 5. Topic-word distributions
  ### ------------------------------------------------------------
  
  psi_out <- init_psi(
    K = K,
    V = V,
    V_set = V_set,
    hyper = hyper
  )
  
  psi      <- psi_out$psi
  psi_star <- psi_out$psi_star
  
  
  ### ------------------------------------------------------------
  ### 6. Correlation parameter
  ### ------------------------------------------------------------
  
  rho <- rbeta(
    1,
    shape1 = hyper$a_rho,
    shape2 = hyper$b_rho
  )
  
  
  ### ------------------------------------------------------------
  ### 7. Assignments + counts
  ### ------------------------------------------------------------
  
  assign <- init_assignments(
    W = W,
    N = N,
    K = K,
    V = V
  )
  
  z    <- assign$z
  n_vk <- assign$n_vk
  n_ik <- assign$n_ik
  
  
  ### ------------------------------------------------------------
  ### 8. Covariance structure
  ### ------------------------------------------------------------
  
  Rstuff <- compute_R(
    rho = rho,
    L = L
  )
  
  R         <- Rstuff$R
  R_inv     <- Rstuff$R_inv
  log_det_R <- Rstuff$log_det_R
  
  
  ### ------------------------------------------------------------
  ### Return full state
  ### ------------------------------------------------------------
  
  return(list(
    
    ## CFA
    mu = mu,
    sigma2 = sigma2,
    phi = phi,
    
    ## Latent factors
    F_mat = F_mat,
    
    ## Topic params
    tau = tau,
    gamma = gamma,
    
    ## Topic-word
    psi = psi,
    psi_star = psi_star,
    
    ## Correlation
    rho = rho,
    
    ## Assignments
    z = z,
    n_vk = n_vk,
    n_ik = n_ik,
    
    ## Covariance
    R = R,
    R_inv = R_inv,
    log_det_R = log_det_R
  ))
}

log_sum_exp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

compute_theta_i <- function(f_i, tau, gamma, lambda_T, B, K, c_val) {
  
  tau_star <- tau
  tau_star[K] <- 0
  
  gamma_star <- gamma
  gamma_star[K] <- 0
  gamma_star[B] <- c_val
  
  logit <- tau_star + gamma_star * f_i[lambda_T]
  
  return(logit - log_sum_exp(logit))
}

log_p_f_i <- function(
    f_i,
    i,
    Y,
    mu, phi, sigma2,
    lambda_S,
    z_i,
    tau_star, gamma_star,
    lambda_T,
    B, c_val,
    R_inv,
    A
) {
  
  J <- ncol(Y)
  K <- length(tau_star)
  M_i <- length(z_i)
  
  
  ### ------------------------------------------------------------
  ### 1. phi*
  ### ------------------------------------------------------------
  
  phi_star <- phi
  phi_star[A] <- 1
  
  
  ### ------------------------------------------------------------
  ### 2. PRIOR: -1/2 f' R^{-1} f
  ### ------------------------------------------------------------
  
  prior_term <- -0.5 * sum(f_i * (R_inv %*% f_i))
  
  
  ### ------------------------------------------------------------
  ### 3. CFA TERM
  ### ------------------------------------------------------------
  
  l_vec <- lambda_S
  
  f_expanded <- f_i[l_vec]
  
  resid <- Y[i, ] - mu - phi_star * f_expanded
  
  cfa_term <- -0.5 * sum((resid^2) / sigma2)
  
  ### ------------------------------------------------------------
  ### 4. LINEAR TERM: sum_m gamma* f
  ### ------------------------------------------------------------
  
  linear_term <- 0
  
  for (m in seq_len(M_i)) {
    k <- z_i[m]
    l_k <- lambda_T[k]
    
    linear_term <- linear_term + gamma_star[k] * f_i[l_k]
  }
  
  
  ### ------------------------------------------------------------
  ### 5. LOG-SUM-EXP TERM
  ### ------------------------------------------------------------
  
  logit <- numeric(K)
  
  for (k in seq_len(K)) {
    l_k <- lambda_T[k]
    logit[k] <- tau_star[k] + gamma_star[k] * f_i[l_k]
  }
  
  log_norm <- log_sum_exp(logit)
  
  norm_term <- - M_i * log_norm
  
  
  ### ------------------------------------------------------------
  ### TOTAL
  ### ------------------------------------------------------------
  
  return(prior_term + cfa_term + linear_term + norm_term)
}

rm_update <- function(log_s, alpha, target, kappa_t) {
  log_s + kappa_t * (alpha - target)
}

empirical_cov_update <- function(mean, cov, x, t) {
  
  # mean: previous mean
  # cov: previous covariance (scaled)
  # x: new sample
  # t: iteration index
  
  if (t == 1) {
    return(list(
      mean = x,
      cov = matrix(0, length(x), length(x))
    ))
  }
  
  delta <- x - mean
  mean_new <- mean + delta / t
  
  cov_new <- cov + tcrossprod(delta, x - mean_new)
  
  return(list(
    mean = mean_new,
    cov = cov_new
  ))
}

build_proposal_cov <- function(s, cov, eps = 1e-6) {
  s^2 * cov + eps * diag(nrow(cov))
}

log_p_tau <- function(
    tau,
    F_mat,
    z,
    gamma,
    lambda_T,
    B, c_val,
    hyper
) {
  
  mu_tau     <- hyper$mu_tau
  sigma_tau2 <- hyper$sigma_tau2
  
  N <- nrow(F_mat)
  K <- length(tau)
  
  
  ### ------------------------------------------------------------
  ### Build tau* and gamma*
  ### ------------------------------------------------------------
  
  tau_star <- numeric(K)
  gamma_star <- numeric(K)
  
  for (k in seq_len(K)) {
    
    tau_star[k] <- if (k == K) 0 else tau[k]
    
    if (k == K) {
      gamma_star[k] <- 0
    } else if (k %in% B) {
      gamma_star[k] <- c_val
    } else {
      gamma_star[k] <- gamma[k]
    }
  }
  
  
  ### ------------------------------------------------------------
  ### Likelihood term
  ### ------------------------------------------------------------
  
  loglik <- 0
  
  for (i in seq_len(N)) {
    
    f_i <- F_mat[i, ]
    z_i <- z[[i]]
    M_i <- length(z_i)
    
    ## linear term
    for (m in seq_len(M_i)) {
      k <- z_i[m]
      loglik <- loglik + tau_star[k]
    }
    
    ## log-sum-exp
    logit <- numeric(K)
    
    for (k in seq_len(K)) {
      l_k <- lambda_T[k]
      logit[k] <- tau_star[k] + gamma_star[k] * f_i[l_k]
    }
    
    loglik <- loglik - M_i * log_sum_exp(logit)
  }
  
  
  ### ------------------------------------------------------------
  ### Prior
  ### ------------------------------------------------------------
  
  tau_nb <- tau[1:(K-1)]
  
  prior <- -0.5 * sum((tau_nb - mu_tau)^2) / sigma_tau2
  
  
  return(loglik + prior)
}

log_p_gamma <- function(
    gamma,
    F_mat,
    z,
    tau,
    lambda_T,
    B, c_val,
    hyper
) {
  
  mu_gamma     <- hyper$mu_gamma
  sigma_gamma2 <- hyper$sigma_gamma2
  
  N <- nrow(F_mat)
  K <- length(gamma)
  
  ### ------------------------------------------------------------
  ### Build gamma* and tau*
  ### ------------------------------------------------------------
  
  gamma_star <- numeric(K)
  tau_star   <- numeric(K)
  
  for (k in seq_len(K)) {
    
    tau_star[k] <- if (k == K) 0 else tau[k]
    
    if (k == K) {
      gamma_star[k] <- 0
    } else if (k %in% B) {
      gamma_star[k] <- c_val
    } else {
      gamma_star[k] <- gamma[k]
    }
  }
  
  
  ### ------------------------------------------------------------
  ### Likelihood
  ### ------------------------------------------------------------
  
  loglik <- 0
  
  for (i in seq_len(N)) {
    
    f_i <- F_mat[i, ]
    z_i <- z[[i]]
    M_i <- length(z_i)
    
    ## linear term
    for (m in seq_len(M_i)) {
      k <- z_i[m]
      l_k <- lambda_T[k]
      loglik <- loglik + gamma_star[k] * f_i[l_k]
    }
    
    ## log-sum-exp
    logit <- numeric(K)
    
    for (k in seq_len(K)) {
      l_k <- lambda_T[k]
      logit[k] <- tau_star[k] + gamma_star[k] * f_i[l_k]
    }
    
    loglik <- loglik - M_i * log_sum_exp(logit)
  }
  
  
  ### ------------------------------------------------------------
  ### Prior (only free indices)
  ### ------------------------------------------------------------
  
  free_idx <- setdiff(1:(K-1), B)
  
  prior <- -0.5 * sum((gamma[free_idx] - mu_gamma)^2) / sigma_gamma2
  
  
  return(loglik + prior)
}

#' Simulate Data from the TACFA Model
#'
#' Generates synthetic data from the Topic-Augmented Confirmatory Factor Analysis (TACFA) model,
#' including both continuous outcomes \eqn{Y} and discrete token data \eqn{W}, along with all
#' structural mappings required for model fitting.
#'
#' The function simulates:
#' \itemize{
#'   \item Latent factors \eqn{\boldsymbol{f}_i}
#'   \item Continuous observations \eqn{Y} via a CFA model
#'   \item Topic proportions \eqn{\boldsymbol{\theta}_i}
#'   \item Topic assignments \eqn{z_{im}} and words \eqn{w_{im}}
#' }
#'
#' ## Model
#'
#' Continuous observations:
#'
#' \deqn{
#' y_{ij} \sim \mathcal{N}\left(
#' \mu_j + \phi_j^* f_{i,\lambda_S(j)},
#' \sigma_j^2
#' \right)
#' }
#'
#' Topic proportions:
#'
#' \deqn{
#' \theta_{ik}
#' =
#' \frac{
#' \exp\left\{\tau_k^* + \gamma_k^* f_{i,\lambda_T(k)}\right\}
#' }{
#' \sum_{k'=1}^K
#' \exp\left\{\tau_{k'}^* + \gamma_{k'}^* f_{i,\lambda_T(k')}\right\}
#' }
#' }
#'
#' Topic assignments and words:
#'
#' \deqn{
#' z_{im} \sim \mathrm{Cat}(\boldsymbol{\theta}_i),
#' \quad
#' w_{im} \mid z_{im} \sim \mathrm{Cat}(\boldsymbol{\psi}^*_{z_{im}})
#' }
#'
#' Latent factors:
#'
#' \deqn{
#' \boldsymbol{f}_i \sim \mathcal{N}(\boldsymbol{0}, R(\rho))
#' }
#'
#' where \eqn{R(\rho)} is an equicorrelation matrix:
#'
#' \deqn{
#' R_{ll'} =
#' \begin{cases}
#' 1 & l = l' \\
#' \rho & l \neq l'
#' \end{cases}
#' }
#'
#' Anchor-word structure:
#'
#' \deqn{
#' \psi^*_{vk} =
#' \begin{cases}
#' \eta & v = v_k \\
#' \epsilon & v \in \mathcal{V},\; v \neq v_k \\
#' \psi_{vk} & v \notin \mathcal{V}
#' \end{cases}
#' }
#'
#' ## Generated Structures
#'
#' The function also constructs all indexing objects required for estimation:
#' \itemize{
#'   \item \eqn{\lambda_S}: variable-to-factor mapping
#'   \item \eqn{\lambda_T}: topic-to-factor mapping
#'   \item \eqn{\mathcal{S}_l}, \eqn{\mathcal{T}_l}: partitions
#'   \item \eqn{\mathcal{A}}: anchor CFA items
#'   \item \eqn{\mathcal{B}}: fixed-loading topics
#'   \item \eqn{\mathcal{V}}: anchor words
#' }
#'
#' @param N Number of observations.
#'
#' @param J Number of observed continuous variables.
#'
#' @param K Number of topics.
#'
#' @param L Number of latent factors.
#'
#' @param V Vocabulary size.
#'
#' @param M_min Minimum number of tokens per observation.
#'
#' @param M_max Maximum number of tokens per observation.
#'
#' @param hyper List of hyperparameters containing:
#' \describe{
#'   \item{mu_0, sigma_mu2}{Prior for \eqn{\mu_j}}
#'   \item{mu_phi, sigma_phi2}{Prior for \eqn{\phi_j}}
#'   \item{a_sigma2, b_sigma2}{Inverse-Gamma for \eqn{\sigma_j^2}}
#'   \item{mu_tau, sigma_tau2}{Prior for \eqn{\tau_k}}
#'   \item{mu_gamma, sigma_gamma2}{Prior for \eqn{\gamma_k}}
#'   \item{a_rho, b_rho}{Beta prior for \eqn{\rho}}
#'   \item{beta}{Dirichlet parameters for \eqn{\boldsymbol{\psi}_k}}
#'   \item{eta, epsilon}{Anchor-word probabilities}
#'   \item{c}{Fixed loading value for topics in \eqn{\mathcal{B}}}
#' }
#'
#' @return A list containing:
#' \describe{
#'   \item{Y}{\eqn{N \times J} matrix of continuous observations}
#'   \item{W}{List of token vectors for each observation}
#'   \item{lambda_S}{Mapping from variables to factors}
#'   \item{lambda_T}{Mapping from topics to factors}
#'   \item{A}{Anchor CFA indices}
#'   \item{B}{Topics with fixed \eqn{\gamma_k = c}}
#'   \item{V_set}{Anchor word indices}
#'   \item{S_list}{Partition of variables}
#'   \item{T_list}{Partition of topics}
#'   \item{true}{List of true parameter values:
#'     \itemize{
#'       \item \eqn{\mu}, \eqn{\phi}, \eqn{\sigma^2}
#'       \item \eqn{\tau}, \eqn{\gamma}
#'       \item \eqn{\rho}
#'       \item \eqn{\boldsymbol{\psi}}, \eqn{\boldsymbol{\psi}^*}
#'       \item latent factors \eqn{\boldsymbol{F}}
#'     }
#'   }
#' }
#'
#' @details
#'
#' The simulation proceeds as follows:
#' \enumerate{
#'   \item Construct partitions \eqn{\mathcal{S}_l} and \eqn{\mathcal{T}_l}
#'   \item Sample global parameters from their priors
#'   \item Generate topic-word distributions \eqn{\boldsymbol{\psi}_k}
#'   \item Apply anchor-word transformation to obtain \eqn{\boldsymbol{\psi}_k^*}
#'   \item Sample latent factors \eqn{\boldsymbol{f}_i}
#'   \item Generate continuous data \eqn{Y}
#'   \item Generate token data \eqn{W} via multinomial sampling
#' }
#'
#' The final output can be passed directly into \code{\link{tacfaSM}} for model fitting.
#'
#' @examples
#' \dontrun{
#' sim <- simulate_tacfa_data(
#'   N = 100,
#'   J = 12,
#'   K = 4,
#'   L = 3,
#'   V = 20,
#'   hyper = hyper
#' )
#'
#' fit <- tacfaSM(
#'   Y = sim$Y,
#'   W = sim$W,
#'   lambda_S = sim$lambda_S,
#'   lambda_T = sim$lambda_T,
#'   A = sim$A,
#'   B = sim$B,
#'   V_set = sim$V_set,
#'   S_list = sim$S_list,
#'   T_list = sim$T_list,
#'   K = 4,
#'   L = 3,
#'   V = 20,
#'   hyper = hyper,
#'   mcmc = list(T = 2000, burn = 1000, chains = 1)
#' )
#' }
#'
#' @seealso \code{\link{tacfaSM}}
#'
#' @export
simulate_tacfa_data <- function(
    N = 100,
    J = 12,
    K = 4,
    L = 3,
    V = 20,
    M_min = 5,
    M_max = 15,
    hyper
) {
  
  ### ------------------------------------------------------------
  ### 1. Construct index structures
  ### ------------------------------------------------------------
  
  # Partition J into S_l (each size >= 3)
  S_list <- split(1:J, rep(1:L, length.out = J))
  lambda_S <- rep(1:L, length.out = J)
  
  # Partition K into T_l
  T_list <- split(1:K, rep(1:L, length.out = K))
  lambda_T <- rep(1:L, length.out = K)
  
  # Anchor items A (choose one per factor)
  A <- sapply(S_list, function(s) s[1])
  
  # B: subset of topics with fixed gamma = c
  B <- sample(1:(K-1), size = floor(K/2))
  
  # V_set: one anchor word per topic
  V_set <- sample(1:V, K, replace = FALSE)
  
  ### ------------------------------------------------------------
  ### 2. Sample global parameters
  ### ------------------------------------------------------------
  
  mu <- rnorm(J, hyper$mu_0, sqrt(hyper$sigma_mu2))
  
  phi <- abs(rnorm(J, hyper$mu_phi, sqrt(hyper$sigma_phi2)))
  
  sigma2 <- 1 / rgamma(J, hyper$a_sigma2, hyper$b_sigma2)
  
  tau <- rnorm(K, hyper$mu_tau, sqrt(hyper$sigma_tau2))
  tau[K] <- 0  # baseline
  
  gamma <- rnorm(K, hyper$mu_gamma, sqrt(hyper$sigma_gamma2))
  gamma[B] <- hyper$c
  gamma[K] <- 0
  
  rho <- rbeta(1, hyper$a_rho, hyper$b_rho)
  
  ### Correlation matrix R(ρ)
  R <- matrix(rho, L, L)
  diag(R) <- 1
  
  ### ------------------------------------------------------------
  ### 3. Sample ψ_k (topic-word distributions)
  ### ------------------------------------------------------------
  
  psi <- matrix(0, V, K)
  
  for (k in 1:K) {
    psi[, k] <- as.numeric(rdirichlet(hyper$beta))
  }
  
  # Build ψ*
  psi_star <- matrix(0, V, K)
  
  for (k in 1:K) {
    for (v in 1:V) {
      if (v %in% V_set) {
        vk <- V_set[k]
        psi_star[v, k] <- ifelse(v == vk, hyper$eta, hyper$epsilon)
      } else {
        psi_star[v, k] <- psi[v, k]
      }
    }
    psi_star[, k] <- psi_star[, k] / sum(psi_star[, k])
  }
  
  ### ------------------------------------------------------------
  ### 4. Sample latent factors F_i
  ### ------------------------------------------------------------
  
  F_mat <- MASS::mvrnorm(N, mu = rep(0, L), Sigma = R)
  
  ### ------------------------------------------------------------
  ### 5. Generate Y
  ### ------------------------------------------------------------
  
  Y <- matrix(0, N, J)
  
  for (i in 1:N) {
    for (j in 1:J) {
      
      l <- lambda_S[j]
      
      phi_star <- if (j %in% A) 1 else phi[j]
      
      mean_ij <- mu[j] + phi_star * F_mat[i, l]
      
      Y[i, j] <- rnorm(1, mean_ij, sqrt(sigma2[j]))
    }
  }
  
  ### ------------------------------------------------------------
  ### 6. Generate W (text / categorical observations)
  ### ------------------------------------------------------------
  
  W <- vector("list", N)
  
  for (i in 1:N) {
    
    Mi <- sample(M_min:M_max, 1)
    W[[i]] <- integer(Mi)
    
    # compute θ_i (softmax)
    logit <- numeric(K)
    
    for (k in 1:K) {
      
      if (k == K) {
        tau_star <- 0
        gamma_star <- 0
      } else {
        tau_star <- tau[k]
        gamma_star <- if (k %in% B) hyper$c else gamma[k]
      }
      
      l <- lambda_T[k]
      
      logit[k] <- tau_star + gamma_star * F_mat[i, l]
    }
    
    theta_i <- exp(logit) / sum(exp(logit))
    
    # sample tokens
    for (m in 1:Mi) {
      z_im <- sample_cat(theta_i)
      w_im <- sample_cat(psi_star[, z_im])
      W[[i]][m] <- w_im
    }
  }
  
  ### ------------------------------------------------------------
  ### 7. Return everything needed for tacfaSM
  ### ------------------------------------------------------------
  
  return(list(
    Y = Y,
    W = W,
    lambda_S = lambda_S,
    lambda_T = lambda_T,
    A = A,
    B = B,
    V_set = V_set,
    S_list = S_list,
    T_list = T_list,
    true = list(
      mu = mu,
      phi = phi,
      sigma2 = sigma2,
      tau = tau,
      gamma = gamma,
      rho = rho,
      psi = psi,
      psi_star = psi_star,
      F = F_mat
    )
  ))
}

#' Generate TACFA Data from Fixed Parameters
#'
#' Generates synthetic data from the Topic-Augmented Confirmatory Factor Analysis (TACFA)
#' model using user-specified parameter values.
#'
#' Unlike \code{\link{simulate_tacfa_data}}, this function does not sample parameters
#' from priors. Instead, it treats all parameters as fixed and generates:
#' \itemize{
#'   \item Continuous observations \eqn{Y}
#'   \item Discrete token data \eqn{W}
#'   \item (Optionally) latent factors \eqn{\boldsymbol{F}}
#' }
#'
#' This is useful for controlled experiments, simulation studies, and recovery analyses.
#'
#' ## Model
#'
#' Continuous outcomes:
#'
#' \deqn{
#' y_{ij} \sim \mathcal{N}\left(
#' \mu_j + \phi_j^* f_{i,\lambda_S(j)},
#' \sigma_j^2
#' \right)
#' }
#'
#' Topic proportions:
#'
#' \deqn{
#' \theta_{ik}
#' =
#' \frac{
#' \exp\left\{\tau_k^* + \gamma_k^* f_{i,\lambda_T(k)}\right\}
#' }{
#' \sum_{k'=1}^K
#' \exp\left\{\tau_{k'}^* + \gamma_{k'}^* f_{i,\lambda_T(k')}\right\}
#' }
#' }
#'
#' Topic assignments and tokens:
#'
#' \deqn{
#' z_{im} \sim \mathrm{Cat}(\boldsymbol{\theta}_i),
#' \quad
#' w_{im} \mid z_{im} \sim \mathrm{Cat}(\boldsymbol{\psi}^*_{z_{im}})
#' }
#'
#' Latent factors:
#'
#' \deqn{
#' \boldsymbol{f}_i \sim \mathcal{N}(\boldsymbol{0}, R(\rho))
#' }
#'
#' with equicorrelation covariance:
#'
#' \deqn{
#' R_{ll'} =
#' \begin{cases}
#' 1 & l = l' \\
#' \rho & l \neq l'
#' \end{cases}
#' }
#'
#' Anchor-word transformation:
#'
#' \deqn{
#' \psi^*_{vk} =
#' \begin{cases}
#' \eta & v = v_k \\
#' \epsilon & v \in \mathcal{V},\; v \neq v_k \\
#' \psi_{vk} & v \notin \mathcal{V}
#' \end{cases}
#' }
#'
#' ## Key Features
#'
#' \itemize{
#'   \item Deterministic parameter input (no prior sampling)
#'   \item Optional user-supplied latent factors \eqn{\boldsymbol{F}}
#'   \item Fully compatible output with \code{\link{tacfaSM}}
#'   \item Anchor-word and structural constraints respected
#' }
#'
#' @param N Number of observations.
#'
#' @param J Number of continuous variables.
#'
#' @param K Number of topics.
#'
#' @param L Number of latent factors.
#'
#' @param V Vocabulary size.
#'
#' @param mu Numeric vector of length \eqn{J} (item intercepts).
#'
#' @param phi Numeric vector of length \eqn{J} (factor loadings).
#'
#' @param sigma2 Numeric vector of length \eqn{J} (error variances).
#'
#' @param tau Numeric vector of length \eqn{K} (topic intercepts).
#'
#' @param gamma Numeric vector of length \eqn{K} (topic loadings).
#'
#' @param rho Scalar in \eqn{(0,1)} controlling factor correlation.
#'
#' @param psi Numeric matrix of dimension \eqn{V \times K} representing
#' topic-word distributions before anchor transformation.
#'
#' @param F_mat Optional \eqn{N \times L} matrix of latent factors.
#' If \code{NULL}, factors are sampled from \eqn{\mathcal{N}(0, R(\rho))}.
#'
#' @param lambda_S Integer vector mapping variables to factors.
#'
#' @param lambda_T Integer vector mapping topics to factors.
#'
#' @param A Integer vector of CFA anchor indices (fixes \eqn{\phi_j = 1}).
#'
#' @param B Integer vector of topics with fixed loading \eqn{\gamma_k = c}.
#'
#' @param V_set Integer vector of anchor words (length \eqn{K}).
#'
#' @param M_min Minimum number of tokens per observation.
#'
#' @param M_max Maximum number of tokens per observation.
#'
#' @param eta Probability assigned to anchor word \eqn{v_k}.
#'
#' @param epsilon Small probability for non-anchor words in \eqn{\mathcal{V}}.
#'
#' @param c Fixed loading value for topics in \eqn{\mathcal{B}}.
#'
#' @return A list containing:
#' \describe{
#'   \item{Y}{\eqn{N \times J} matrix of continuous observations}
#'   \item{W}{List of token vectors}
#'   \item{F}{Latent factor matrix}
#'   \item{psi_star}{Transformed topic-word distributions}
#'   \item{R}{Factor covariance matrix}
#' }
#'
#' @details
#'
#' The generation proceeds as follows:
#' \enumerate{
#'   \item Construct covariance matrix \eqn{R(\rho)}
#'   \item Generate latent factors if not provided
#'   \item Apply anchor-word transformation to \eqn{\boldsymbol{\psi}}
#'   \item Generate continuous outcomes \eqn{Y}
#'   \item Generate token data \eqn{W} using multinomial sampling
#' }
#'
#' The softmax computation for \eqn{\boldsymbol{\theta}_i} is stabilized by subtracting
#' \eqn{\max_k} from logits before exponentiation.
#'
#' @examples
#' \dontrun{
#' sim <- generate_tacfa_from_params(
#'   N = 100,
#'   J = 12,
#'   K = 4,
#'   L = 3,
#'   V = 20,
#'   mu = mu,
#'   phi = phi,
#'   sigma2 = sigma2,
#'   tau = tau,
#'   gamma = gamma,
#'   rho = 0.3,
#'   psi = psi,
#'   lambda_S = lambda_S,
#'   lambda_T = lambda_T,
#'   A = A,
#'   B = B,
#'   V_set = V_set,
#'   eta = 0.8,
#'   epsilon = 0.01,
#'   c = 1
#' )
#' }
#'
#' @seealso \code{\link{simulate_tacfa_data}}, \code{\link{tacfaSM}}
#'
#' @export
generate_tacfa_from_params <- function(
    N,
    J, K, L, V,
    
    # data-generating parameters
    mu,
    phi,
    sigma2,
    tau,
    gamma,
    rho,
    psi,
    
    # latent factors (optional: if NULL → generated from R)
    F_mat = NULL,
    
    # structure
    lambda_S,
    lambda_T,
    A,
    B,
    V_set,
    
    # text length
    M_min = 5,
    M_max = 15,
    
    # hyper constants used in ψ*
    eta,
    epsilon,
    c
) {
  
  ### ------------------------------------------------------------
  ### 0. Basic checks
  ### ------------------------------------------------------------
  
  if (length(mu) != J) stop("mu length mismatch")
  if (length(phi) != J) stop("phi length mismatch")
  if (length(sigma2) != J) stop("sigma2 length mismatch")
  
  if (length(tau) != K) stop("tau length mismatch")
  if (length(gamma) != K) stop("gamma length mismatch")
  
  if (!all(V_set %in% 1:V)) stop("V_set invalid")
  if (length(V_set) != K) stop("|V_set| must equal K")
  
  ### ------------------------------------------------------------
  ### 1. Build R(ρ)
  ### ------------------------------------------------------------
  
  R <- matrix(rho, L, L)
  diag(R) <- 1
  
  ### ------------------------------------------------------------
  ### 2. Generate latent factors if not supplied
  ### ------------------------------------------------------------
  
  if (is.null(F_mat)) {
    F_mat <- MASS::mvrnorm(N, mu = rep(0, L), Sigma = R)
  }
  
  ### ------------------------------------------------------------
  ### 3. Build ψ*
  ### ------------------------------------------------------------
  
  psi_star <- matrix(0, V, K)
  
  for (k in 1:K) {
    for (v in 1:V) {
      
      if (v %in% V_set) {
        vk <- V_set[k]
        psi_star[v, k] <- ifelse(v == vk, eta, epsilon)
      } else {
        psi_star[v, k] <- psi[v, k]
      }
    }
    psi_star[, k] <- psi_star[, k] / sum(psi_star[, k])
  }
  
  ### ------------------------------------------------------------
  ### 4. Generate Y
  ### ------------------------------------------------------------
  
  Y <- matrix(0, N, J)
  
  for (i in 1:N) {
    for (j in 1:J) {
      
      l <- lambda_S[j]
      
      phi_star <- if (j %in% A) 1 else phi[j]
      
      mean_ij <- mu[j] + phi_star * F_mat[i, l]
      
      Y[i, j] <- rnorm(1, mean_ij, sqrt(sigma2[j]))
    }
  }
  
  ### ------------------------------------------------------------
  ### 5. Generate W
  ### ------------------------------------------------------------
  
  W <- vector("list", N)
  
  for (i in 1:N) {
    
    Mi <- sample(M_min:M_max, 1)
    W[[i]] <- integer(Mi)
    
    ### ---- compute θ_i ----
    logit <- numeric(K)
    
    for (k in 1:K) {
      
      if (k == K) {
        tau_star <- 0
        gamma_star <- 0
      } else {
        tau_star <- tau[k]
        gamma_star <- if (k %in% B) c else gamma[k]
      }
      
      l <- lambda_T[k]
      
      logit[k] <- tau_star + gamma_star * F_mat[i, l]
    }
    
    # numerical stability
    logit <- logit - max(logit)
    
    theta_i <- exp(logit) / sum(exp(logit))
    
    ### ---- sample tokens ----
    for (m in 1:Mi) {
      
      z_im <- sample_cat(theta_i)
      w_im <- sample_cat(psi_star[, z_im])
      
      W[[i]][m] <- w_im
    }
  }
  
  ### ------------------------------------------------------------
  ### 6. Return
  ### ------------------------------------------------------------
  
  return(list(
    Y = Y,
    W = W,
    F = F_mat,
    psi_star = psi_star,
    R = R
  ))
}

log_p_rho <- function(rho, F_mat, hyper, L) {
  
  ### ------------------------------------------------------------
  ### 0. Domain check
  ### ------------------------------------------------------------
  
  if (rho <= 0 || rho >= 1) return(-Inf)
  
  N <- nrow(F_mat)
  
  ### ------------------------------------------------------------
  ### 1. Build R(ρ)
  ### ------------------------------------------------------------
  
  # Compound symmetry matrix
  R <- matrix(rho, L, L)
  diag(R) <- 1
  
  ### ------------------------------------------------------------
  ### 2. Check positive definiteness
  ### ------------------------------------------------------------
  
  eigvals <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigvals <= 0)) return(-Inf)
  
  ### ------------------------------------------------------------
  ### 3. Inverse + log determinant
  ### ------------------------------------------------------------
  
  R_inv <- solve(R)
  log_det_R <- as.numeric(determinant(R, logarithm = TRUE)$modulus)
  
  ### ------------------------------------------------------------
  ### 4. Likelihood: F_i ~ N(0, R)
  ### ------------------------------------------------------------
  
  quad_terms <- rowSums((F_mat %*% R_inv) * F_mat)
  
  loglik <- -0.5 * N * log_det_R - 0.5 * sum(quad_terms)
  
  ### ------------------------------------------------------------
  ### 5. Prior: Beta(a_rho, b_rho)
  ### ------------------------------------------------------------
  
  logprior <- dbeta(rho, hyper$a_rho, hyper$b_rho, log = TRUE)
  
  ### ------------------------------------------------------------
  ### 6. Return log posterior (up to constant)
  ### ------------------------------------------------------------
  
  return(loglik + logprior)
}

#' Construct Hyperparameters for the TACFA Model
#'
#' Creates a validated list of hyperparameters for the
#' Topic-Augmented Confirmatory Factor Analysis (TACFA) model.
#'
#' This function provides a convenient and consistent way to specify priors
#' used in both simulation and estimation functions such as
#' \code{\link{simulate_tacfa_data}} and \code{\link{tacfaSM}}.
#'
#' ## Priors
#'
#' The following priors are defined:
#'
#' \strong{CFA component}
#'
#' \deqn{
#' \mu_j \sim \mathcal{N}(\mu_0, \sigma_\mu^2)
#' }
#'
#' \deqn{
#' \phi_j \sim \mathcal{N}^+(\mu_\phi, \sigma_\phi^2)
#' }
#'
#' \deqn{
#' \sigma_j^2 \sim \mathrm{InvGamma}(a_{\sigma^2}, b_{\sigma^2})
#' }
#'
#' \strong{Topic regression component}
#'
#' \deqn{
#' \tau_k \sim \mathcal{N}(\mu_\tau, \sigma_\tau^2)
#' }
#'
#' \deqn{
#' \gamma_k \sim \mathcal{N}(\mu_\gamma, \sigma_\gamma^2)
#' }
#'
#' \strong{Latent factor correlation}
#'
#' \deqn{
#' \rho \sim \mathrm{Beta}(a_\rho, b_\rho)
#' }
#'
#' \strong{Topic-word distributions}
#'
#' \deqn{
#' \boldsymbol{\psi}_k \sim \mathrm{Dirichlet}(\boldsymbol{\beta})
#' }
#'
#' \strong{Anchor-word structure}
#'
#' \deqn{
#' \psi^*_{vk} =
#' \begin{cases}
#' \eta & v = v_k \\
#' \epsilon & v \in \mathcal{V},\; v \neq v_k \\
#' (1 - \eta - (K-1)\epsilon)\,\tilde{\psi}_{vk} & v \notin \mathcal{V}
#' \end{cases}
#' }
#'
#' where \eqn{\eta \gg \epsilon} controls anchor strength.
#'
#' @param mu_0 Mean of the prior for \eqn{\mu_j}.
#'
#' @param sigma_mu2 Variance of the prior for \eqn{\mu_j}.
#'
#' @param mu_phi Mean of the prior for \eqn{\phi_j}.
#'
#' @param sigma_phi2 Variance of the prior for \eqn{\phi_j}.
#'
#' @param a_sigma2 Shape parameter of the inverse-gamma prior for \eqn{\sigma_j^2}.
#'
#' @param b_sigma2 Scale parameter of the inverse-gamma prior for \eqn{\sigma_j^2}.
#'
#' @param mu_tau Mean of the prior for \eqn{\tau_k}.
#'
#' @param sigma_tau2 Variance of the prior for \eqn{\tau_k}.
#'
#' @param mu_gamma Mean of the prior for \eqn{\gamma_k}.
#'
#' @param sigma_gamma2 Variance of the prior for \eqn{\gamma_k}.
#'
#' @param a_rho Shape parameter of the Beta prior for \eqn{\rho}.
#'
#' @param b_rho Shape parameter of the Beta prior for \eqn{\rho}.
#'
#' @param beta Numeric vector of Dirichlet parameters for topic-word distributions.
#' Must be strictly positive.
#'
#' @param eta Probability mass assigned to anchor words.
#'
#' @param epsilon Small probability for non-anchor words within the anchor set.
#'
#' @param c Fixed loading value used for topics in \eqn{\mathcal{B}}.
#'
#' @param V Vocabulary size. Only used if \code{beta = NULL}, in which case
#' \code{beta} is initialized as a symmetric Dirichlet vector.
#'
#' @return A named list containing all hyperparameters required by the TACFA model:
#' \describe{
#'   \item{mu_0, sigma_mu2}{CFA intercept prior}
#'   \item{mu_phi, sigma_phi2}{CFA loading prior}
#'   \item{a_sigma2, b_sigma2}{Variance prior}
#'   \item{mu_tau, sigma_tau2}{Topic intercept prior}
#'   \item{mu_gamma, sigma_gamma2}{Topic loading prior}
#'   \item{a_rho, b_rho}{Correlation prior}
#'   \item{beta}{Dirichlet parameters}
#'   \item{eta, epsilon}{Anchor-word parameters}
#'   \item{c}{Fixed topic loading}
#' }
#'
#' @details
#'
#' If \code{beta} is not provided, it is initialized as:
#'
#' \deqn{
#' \beta_v = 0.1 \quad \forall v \in \{1,\dots,V\}
#' }
#'
#' This corresponds to a symmetric Dirichlet prior encouraging moderately sparse
#' topic distributions.
#'
#' Note that full validity of the anchor-word construction requires:
#'
#' \deqn{
#' \eta + (K - 1)\epsilon < 1
#' }
#'
#' which depends on \eqn{K} and is therefore checked during model execution.
#'
#' @examples
#' hyper <- make_tacfa_hyper(
#'   mu_0 = 0,
#'   sigma_mu2 = 1,
#'   beta = rep(0.1, 20),
#'   eta = 0.7,
#'   epsilon = 0.01
#' )
#'
#' @seealso \code{\link{tacfaSM}}, \code{\link{simulate_tacfa_data}}
#'
#' @export
make_tacfa_hyper <- function(
    # CFA priors
  mu_0 = 0,
  sigma_mu2 = 1,
  mu_phi = 1,
  sigma_phi2 = 0.5,
  
  # Variance prior
  a_sigma2 = 2,
  b_sigma2 = 2,
  
  # Topic regression priors
  mu_tau = 0,
  sigma_tau2 = 1,
  mu_gamma = 0,
  sigma_gamma2 = 1,
  
  # Correlation prior
  a_rho = 2,
  b_rho = 2,
  
  # Dirichlet prior
  beta = NULL,   # must be length V (checked below)
  
  # Anchor structure
  eta = 0.6,
  epsilon = 0.01,
  
  # Fixed gamma value
  c = 1,
  
  # Vocabulary size (only needed if beta is NULL)
  V = NULL
) {
  
  ## ------------------------------------------------------------
  ## Validate / build beta
  ## ------------------------------------------------------------
  
  if (is.null(beta)) {
    if (is.null(V)) {
      stop("Either beta or V must be provided")
    }
    beta <- rep(0.1, V)
  }
  
  if (!is.numeric(beta) || any(beta <= 0)) {
    stop("beta must be a positive numeric vector")
  }
  
  V <- length(beta)
  
  ## ------------------------------------------------------------
  ## Validate anchor mass
  ## ------------------------------------------------------------
  
  remaining_mass <- 1 - eta - (length(beta) > 0) * 0  # just structural
  
  if (eta <= 0 || epsilon <= 0) {
    stop("eta and epsilon must be positive")
  }
  
  ## NOTE: full check requires K, so defer strict check to model
  
  ## ------------------------------------------------------------
  ## Build list
  ## ------------------------------------------------------------
  
  hyper <- list(
    mu_0 = mu_0,
    sigma_mu2 = sigma_mu2,
    
    mu_phi = mu_phi,
    sigma_phi2 = sigma_phi2,
    
    a_sigma2 = a_sigma2,
    b_sigma2 = b_sigma2,
    
    mu_tau = mu_tau,
    sigma_tau2 = sigma_tau2,
    
    mu_gamma = mu_gamma,
    sigma_gamma2 = sigma_gamma2,
    
    a_rho = a_rho,
    b_rho = b_rho,
    
    beta = beta,
    
    eta = eta,
    epsilon = epsilon,
    
    c = c
  )
  
  return(hyper)
}
