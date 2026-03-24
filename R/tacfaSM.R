#' @useDynLib tacfa, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom posterior as_draws_array
#' @importFrom MASS mvrnorm
#' @importFrom stats dbeta rbeta rgamma rnorm runif
#' @importFrom utils modifyList
#' @importFrom abind abind
NULL

#' Topic-Augmented Confirmatory Factor Analysis with Structured Multinomial Component
#'
#' Fits the joint CFA–topic model using a Gibbs sampler with adaptive
#' Metropolis–Hastings steps for continuous latent parameters.
#'
#' This model combines:
#' \itemize{
#'   \item A Gaussian confirmatory factor analysis (CFA) model for continuous data \eqn{Y}
#'   \item A multinomial topic model (LDA-style) for discrete observations \eqn{W}
#'   \item Shared latent factors \eqn{\boldsymbol{f}_i} linking both components
#' }
#'
#' ## Model
#'
#' The model is defined as:
#'
#' \deqn{
#' y_{ij} \mid \mu_j, \phi_j, \boldsymbol{f}_i, \sigma_j^2
#' \sim \mathcal{N}\left(
#' \mu_j + \phi_j^* f_{i,\lambda_S(j)},
#' \sigma_j^2
#' \right)
#' }
#'
#' \deqn{
#' z_{im} \sim \mathrm{Cat}(\boldsymbol{\theta}_i),
#' \quad
#' w_{im} \mid z_{im} \sim \mathrm{Cat}(\boldsymbol{\psi}^*_{z_{im}})
#' }
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
#' Latent factors follow:
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
#' ## Inference
#'
#' Posterior inference is performed via a Gibbs sampler with:
#' \itemize{
#'   \item Closed-form updates for \eqn{\mu_j}, \eqn{\sigma_j^2}, \eqn{\phi_j}, and \eqn{\psi_k}
#'   \item Categorical sampling for topic assignments \eqn{z_{im}}
#'   \item Adaptive Metropolis–Hastings updates for:
#'   \itemize{
#'     \item latent factors \eqn{\boldsymbol{f}_i}
#'     \item topic intercepts \eqn{\boldsymbol{\tau}}
#'     \item topic loadings \eqn{\boldsymbol{\gamma}}
#'     \item correlation parameter \eqn{\rho}
#'   }
#' }
#'
#' Adaptive proposals use Robbins–Monro scaling during burn-in.
#'
#' @param Y Numeric matrix of dimension \eqn{N \times J} containing continuous observations.
#'
#' @param W List of length \eqn{N}. Each element \eqn{W[[i]]} is a vector of
#' word indices in \eqn{\{1,\dots,V\}}.
#'
#' @param lambda_S Integer vector of length \eqn{J}. Maps each observed variable
#' to a latent factor: \eqn{\lambda_S(j) \in \{1,\dots,L\}}.
#'
#' @param lambda_T Integer vector of length \eqn{K}. Maps each topic to a latent factor:
#' \eqn{\lambda_T(k)}.
#'
#' @param A Integer vector specifying indices of items with fixed loading \eqn{\phi_j = 1}.
#'
#' @param B Integer vector specifying topics with fixed loading \eqn{\gamma_k = c}.
#'
#' @param V_set Integer vector of length \eqn{K} specifying anchor words.
#' Must be unique and satisfy \eqn{|V_{set}| = K}.
#'
#' @param S_list List of length \eqn{L}. Each element contains indices of observed
#' variables assigned to factor \eqn{l}.
#'
#' @param T_list List of length \eqn{L}. Each element contains indices of topics
#' assigned to factor \eqn{l}.
#'
#' @param K Number of topics.
#'
#' @param L Number of latent factors.
#'
#' @param V Vocabulary size.
#'
#' @param hyper List of hyperparameters containing:
#' \describe{
#'   \item{mu_0, sigma_mu2}{Prior for \eqn{\mu_j}}
#'   \item{mu_phi, sigma_phi2}{Prior for \eqn{\phi_j}}
#'   \item{a_sigma2, b_sigma2}{Inverse-Gamma prior for \eqn{\sigma_j^2}}
#'   \item{mu_tau, sigma_tau2}{Prior for \eqn{\tau_k}}
#'   \item{mu_gamma, sigma_gamma2}{Prior for \eqn{\gamma_k}}
#'   \item{a_rho, b_rho}{Beta prior for \eqn{\rho}}
#'   \item{beta}{Dirichlet prior for topic-word distributions}
#'   \item{eta, epsilon}{Anchor-word probabilities}
#'   \item{c}{Fixed loading value for \eqn{k \in B}}
#' }
#'
#' @param mcmc List of MCMC settings:
#' \describe{
#'   \item{T}{Total iterations}
#'   \item{burn}{Burn-in iterations}
#'   \item{chains}{Number of chains}
#'   \item{delta_f, delta_tau, delta_gamma}{Initial proposal scales}
#'   \item{sigma_rho}{Proposal sd for transformed \eqn{\rho}}
#'   \item{adapt}{Logical; enable adaptation}
#'   \item{a}{Robbins–Monro exponent}
#' }
#'
#' @param init Optional list of initial values (currently unused).
#'
#' @param control List of control options:
#' \describe{
#'   \item{verbose}{Print progress}
#'   \item{store_z}{Store latent assignments}
#'   \item{check_inputs}{Validate inputs}
#'   \item{return_draws}{Return posterior draws}
#' }
#'
#' @return An object of class \code{"tacfaSM"} containing:
#' \describe{
#'   \item{draws}{Posterior draws (if requested)}
#'   \item{raw}{Arrays of stored samples}
#'   \item{dims}{Model dimensions}
#'   \item{call}{Function call}
#' }
#'
#' @details
#'
#' The sampler proceeds in the following steps each iteration:
#' \enumerate{
#'   \item Sample topic assignments \eqn{z_{im}}
#'   \item Sample topic-word distributions \eqn{\boldsymbol{\psi}_k}
#'   \item Sample \eqn{\mu_j}, \eqn{\sigma_j^2}, \eqn{\phi_j}
#'   \item Sample latent factors \eqn{\boldsymbol{f}_i} via adaptive MH
#'   \item Sample \eqn{\boldsymbol{\tau}} and \eqn{\boldsymbol{\gamma}} via adaptive MH
#'   \item Sample \eqn{\rho} using a logit transform
#' }
#'
#' Adaptation is performed during burn-in and then frozen to ensure ergodicity.
#'
#' @examples
#' \dontrun{
#' fit <- tacfaSM(
#'   Y = Y,
#'   W = W,
#'   lambda_S = lambda_S,
#'   lambda_T = lambda_T,
#'   A = A,
#'   B = B,
#'   V_set = V_set,
#'   S_list = S_list,
#'   T_list = T_list,
#'   K = K,
#'   L = L,
#'   V = V,
#'   hyper = hyper,
#'   mcmc = list(T = 2000, burn = 1000, chains = 2)
#' )
#' }
#'
#' @references
#' This model combines ideas from confirmatory factor analysis and latent Dirichlet
#' allocation with shared latent structure.
#'
#'
#' @export
tacfaSM <- function(
    Y, W,
    lambda_S, lambda_T,
    A, B, V_set,
    S_list, T_list,
    K, L, V,
    hyper,
    mcmc,
    init = NULL,
    control = list()
) {
  
  ### ------------------------------------------------------------
  ### 0. Control defaults
  ### ------------------------------------------------------------
  
  control_default <- list(
    verbose = TRUE,
    store_z = FALSE,
    check_inputs = TRUE,
    return_draws = TRUE
  )
  
  control <- modifyList(control_default, control)
  
  
  ### ------------------------------------------------------------
  ### 1. Basic data checks
  ### ------------------------------------------------------------
  
  if (!is.matrix(Y)) stop("Y must be a matrix")
  if (!is.list(W)) stop("W must be a list")
  
  N <- nrow(Y)
  J <- ncol(Y)
  
  if (length(W) != N) stop("W must have length N")
  
  M <- sapply(W, length)
  
  
  ### ------------------------------------------------------------
  ### 2. Input validation
  ### ------------------------------------------------------------
  
  if (control$check_inputs) {
    
    if (length(lambda_S) != J)
      stop("lambda_S must have length J")
    if (!all(lambda_S %in% seq_len(L)))
      stop("lambda_S invalid")
    
    if (length(lambda_T) != K)
      stop("lambda_T must have length K")
    if (!all(lambda_T %in% seq_len(L)))
      stop("lambda_T invalid")
    
    if (!all(A %in% seq_len(J)))
      stop("A invalid")
    
    if (!all(B %in% seq_len(K)))
      stop("B invalid")
    
    if (!all(V_set %in% seq_len(V)))
      stop("V_set invalid")
    if (length(V_set) != K)
      stop("|V_set| must equal K")
    if (length(unique(V_set)) != K)
      stop("V_set must be unique")
    
    if (!identical(unname(sort(unlist(S_list))), seq_len(J)))
      stop("S_list invalid")
    
    if (!identical(unname(sort(unlist(T_list))), seq_len(K)))
      stop("T_list invalid")
    
    for (i in seq_len(N)) {
      if (!all(W[[i]] %in% seq_len(V))) {
        stop(sprintf("W[[%d]] invalid", i))
      }
    }
  }
  
  
  ### ------------------------------------------------------------
  ### 3. Hyperparameter validation
  ### ------------------------------------------------------------
  
  required_hyper <- c(
    "mu_0","sigma_mu2",
    "mu_phi","sigma_phi2",
    "a_sigma2","b_sigma2",
    "mu_tau","sigma_tau2",
    "mu_gamma","sigma_gamma2",
    "a_rho","b_rho",
    "beta","eta","epsilon","c"
  )
  
  missing_hyper <- setdiff(required_hyper, names(hyper))
  if (length(missing_hyper) > 0) {
    stop(paste("Missing hyperparameters:", paste(missing_hyper, collapse = ", ")))
  }
  
  
  ### ------------------------------------------------------------
  ### 4. MCMC settings
  ### ------------------------------------------------------------
  
  T_iter <- mcmc$T
  burn   <- mcmc$burn
  chains <- mcmc$chains
  
  if (burn >= T_iter) stop("burn must be < T")
  
  n_store <- T_iter - burn
  
  
  ### ------------------------------------------------------------
  ### 5. Storage
  ### ------------------------------------------------------------
  
  mu_store     <- array(NA_real_, c(n_store, chains, J))
  sigma2_store <- array(NA_real_, c(n_store, chains, J))
  phi_store    <- array(NA_real_, c(n_store, chains, J))
  
  F_store      <- array(NA_real_, c(n_store, chains, N, L))
  
  tau_store    <- array(NA_real_, c(n_store, chains, K))
  gamma_store  <- array(NA_real_, c(n_store, chains, K))
  
  rho_store    <- array(NA_real_, c(n_store, chains))
  
  
  ### ------------------------------------------------------------
  ### 6. Run chains
  ### ------------------------------------------------------------
  
  for (ch in seq_len(chains)) {
    
    if (control$verbose) cat("Chain", ch, "\n")
    
    
    ### ---------------- INIT ----------------
    
    state <- initialize_parameters(
      Y, W,
      lambda_S, lambda_T,
      A, B, V_set,
      S_list, T_list,
      K, L, V,
      hyper
    )
    
    mu     <- state$mu
    sigma2 <- state$sigma2
    phi    <- state$phi
    
    F_mat  <- state$F_mat
    
    tau    <- state$tau
    gamma  <- state$gamma
    
    psi      <- state$psi
    psi_star <- state$psi_star
    
    rho    <- state$rho
    
    z      <- state$z
    n_vk   <- state$n_vk
    n_ik   <- state$n_ik
    
    R        <- state$R
    R_inv    <- state$R_inv
    log_detR <- state$log_det_R
    
    
    ### ----------- ADAPTIVE STATE -----------
    
    adapt_f <- lapply(seq_len(N), function(i) {
      list(
        log_s = log(mcmc$delta_f),
        mean  = F_mat[i, ],
        cov   = matrix(0, L, L),
        t     = 1
      )
    })
    
    adapt_tau <- list(
      log_s = log(mcmc$delta_tau),
      mean  = tau[1:(K-1)],
      cov   = matrix(0, K-1, K-1),
      t     = 1
    )
    
    free_gamma_idx <- setdiff(1:(K-1), B)
    
    adapt_gamma <- list(
      log_s = log(mcmc$delta_gamma),
      mean  = gamma[free_gamma_idx],
      cov   = matrix(0, length(free_gamma_idx), length(free_gamma_idx)),
      t     = 1
    )
    
    adapt_rho <- list(
      log_s = log(mcmc$sigma_rho)
    )
    
    
    store_idx <- 1
    
    
    ### ---------------- MCMC ----------------
    
    for (t in seq_len(T_iter)) {
      
      ## 1 z
      res <- sample_z(W,z,F_mat,tau,gamma,psi_star,lambda_T,B,hyper$c,n_vk,n_ik)
      z<-res$z; n_vk<-res$n_vk; n_ik<-res$n_ik
      
      ## 2 psi
      psi_out <- sample_psi(n_vk,V_set,hyper)
      psi<-psi_out$psi; psi_star<-psi_out$psi_star
      
      ## 3 mu
      mu <- sample_mu(Y,F_mat,phi,sigma2,lambda_S,A,hyper)
      
      ## 4 sigma2
      sigma2 <- sample_sigma2(Y,F_mat,mu,phi,lambda_S,A,hyper)
      
      ## 5 phi
      phi <- sample_phi(Y,F_mat,mu,sigma2,lambda_S,A,hyper)
      
      ## 6 F
      resF <- sample_F(F_mat,adapt_f,t,mcmc$adapt,mcmc$a,
                       Y,mu,phi,sigma2,lambda_S,z,
                       tau,gamma,lambda_T,B,hyper$c,
                       R_inv,A)
      F_mat<-resF$F_mat; adapt_f<-resF$adapt_f
      
      ## 7 tau
      resT <- sample_tau(tau,adapt_tau,t,mcmc$adapt,mcmc$a,
                         F_mat,z,gamma,lambda_T,B,hyper$c,hyper)
      tau<-resT$tau; adapt_tau<-resT$adapt_tau
      
      ## 8 gamma
      resG <- sample_gamma(gamma,adapt_gamma,t,mcmc$adapt,mcmc$a,
                           F_mat,z,tau,lambda_T,B,hyper$c,hyper)
      gamma<-resG$gamma; adapt_gamma<-resG$adapt_gamma
      
      ## 9 rho
      resR <- sample_rho(rho,adapt_rho,t,mcmc$adapt,mcmc$a,
                         F_mat,hyper,L)
      rho<-resR$rho
      adapt_rho<-resR$adapt_rho
      R<-resR$R; R_inv<-resR$R_inv; log_detR<-resR$log_det_R
      
      
      ## store
      if (t > burn) {
        mu_store[store_idx,ch,]<-mu
        sigma2_store[store_idx,ch,]<-sigma2
        phi_store[store_idx,ch,]<-phi
        F_store[store_idx,ch,,]<-F_mat
        tau_store[store_idx,ch,]<-tau
        gamma_store[store_idx,ch,]<-gamma
        rho_store[store_idx,ch]<-rho
        store_idx<-store_idx+1
      }
    }
  }
  
  
  ### ------------------------------------------------------------
  ### 7. Posterior draws
  ### ------------------------------------------------------------
  
  draws <- if (control$return_draws) {
    build_draws_chains(
      mu_store,sigma2_store,phi_store,
      F_store,tau_store,gamma_store,rho_store
    )
  } else NULL
  
  
  ### ------------------------------------------------------------
  ### 8. Return
  ### ------------------------------------------------------------
  
  out <- list(
    draws = draws,
    raw = list(
      mu = mu_store,
      sigma2 = sigma2_store,
      phi = phi_store,
      F = F_store,
      tau = tau_store,
      gamma = gamma_store,
      rho = rho_store
    ),
    dims = list(N=N,J=J,K=K,L=L,V=V),
    call = match.call()
  )
  
  class(out) <- "tacfaSM"
  
  return(out)
}