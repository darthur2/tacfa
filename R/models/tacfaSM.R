library(R6)

TacfaSM <- R6Class(
  "TacfaSM",
  
  public = list(
    
    # -----------------------------
    # Dimensions
    # -----------------------------
    L = NULL,
    K = NULL,
    V = NULL,
    
    # -----------------------------
    # Seed structure
    # -----------------------------
    seed_index = NULL, # vector
    psi_seed = NULL,   # matrix
    V_free = NULL,     # vector
    
    # -----------------------------
    # Hyperparameters
    # -----------------------------
    mu0 = NULL,
    sigma_mu2 = NULL,
    
    mu_phi = NULL,
    sigma_phi2 = NULL,
    
    a_sigma2 = NULL,
    b_sigma2 = NULL,
    
    a_rho = NULL,
    b_rho = NULL,
    
    mu_delta = NULL,
    sigma_delta2 = NULL,
    
    mu_gamma = NULL,
    sigma_gamma2 = NULL,
    
    beta_V = NULL,
    
    # -----------------------------
    # Initialize (NO DATA)
    # -----------------------------
    initialize = function(
    L, K, V,
    seed_index = integer(0),
    psi_seed = NULL,
    beta_V = NULL,
    
    mu0 = 0,
    sigma_mu2 = 1,
    mu_phi = 0,
    sigma_phi2 = 1,
    a_sigma2 = 2,
    b_sigma2 = 2,
    a_rho = 2,
    b_rho = 2,
    mu_delta = 0,
    sigma_delta2 = 1,
    mu_gamma = 0,
    sigma_gamma2 = 1
    ) {
      
      self$L <- L
      self$K <- K
      self$V <- V
      
      self$seed_index <- seed_index
      self$psi_seed <- psi_seed
      self$V_free <- setdiff(1:V, seed_index)
      
      self$beta_V <- beta_V
      
      # Hyperparameters
      self$mu0 <- mu0
      self$sigma_mu2 <- sigma_mu2
      
      self$mu_phi <- mu_phi
      self$sigma_phi2 <- sigma_phi2
      
      self$a_sigma2 <- a_sigma2
      self$b_sigma2 <- b_sigma2
      
      self$a_rho <- a_rho
      self$b_rho <- b_rho
      
      self$mu_delta <- mu_delta
      self$sigma_delta2 <- sigma_delta2
      
      self$mu_gamma <- mu_gamma
      self$sigma_gamma2 <- sigma_gamma2
      
      # Validation
      if (length(seed_index) > 0) {
        if (is.null(psi_seed)) {
          stop("psi_seed must be provided if seed_index is non-empty")
        }
        if (ncol(psi_seed) != K) {
          stop("psi_seed must have K columns")
        }
      }
    },
    
    # -----------------------------
    # Construct ψ
    # -----------------------------
    build_psi = function(psi_free) {
      
      psi <- matrix(0, self$V, self$K)
      
      for (k in 1:self$K) {
        
        c_k <- if (length(self$seed_index) > 0) {
          sum(self$psi_seed[, k])
        } else 0
        
        if (c_k >= 1) {
          stop("Seed probabilities must sum to < 1")
        }
        
        psi[self$V_free, k] <- (1 - c_k) * psi_free[, k]
        
        if (length(self$seed_index) > 0) {
          psi[self$seed_index, k] <- self$psi_seed[, k]
        }
      }
      
      psi
    },
    
    # -----------------------------
    # Softmax
    # -----------------------------
    softmax = function(x) {
      x <- x - max(x)
      exp(x) / sum(exp(x))
    },
    
    # -----------------------------
    # Compute θ_i
    # -----------------------------
    compute_theta = function(f_i, delta, Gamma) {
      eta <- delta + as.vector(f_i %*% Gamma)
      self$softmax(eta)
    },
    
    # -----------------------------
    # Fit (ENGINE-BASED)
    # -----------------------------
    fit = function(Y, W, lambda, engine = NULL, ...) {
      
      # Default engine = Gibbs
      if (is.null(engine)) {
        engine <- TacfaSMGibbs$new(...)
      }
      
      # Attach model + data
      engine$setup(
        model = self,
        Y = Y,
        W = W,
        lambda = lambda
      )
      
      # Run inference
      engine$run()
      
      # Return fitted state
      return(engine$get_samples())
    },
    
    generate_test_data = function(
    N,
    J,
    M = 50,
    seed = NULL
    ) {
      
      if (!is.null(seed)) set.seed(seed)
      
      L <- self$L
      K <- self$K
      V <- self$V
      
      # -----------------------------
      # 1. Generate CFA structure
      # -----------------------------
      
      lambda <- rep(1:L, each = J/L)
      
      mu <- sample(c(-1, 1), J, TRUE)
      phi <- rep(1, J)  # positive loadings
      sigma2 <- rep(1.5, J)
      
      # -----------------------------
      # 2. Generate latent factors
      # -----------------------------
      
      rho <- 0.3
      R <- matrix(rho, L, L)
      diag(R) <- 1
      
      FF <- MASS::mvrnorm(N, mu = rep(0, L), Sigma = R)
      
      # -----------------------------
      # 3. Generate Y
      # -----------------------------
      
      Y <- matrix(0, N, J)
      
      for (i in 1:N) {
        for (j in 1:J) {
          mean_ij <- mu[j] + phi[j] * FF[i, lambda[j]]
          Y[i, j] <- rnorm(1, mean_ij, sqrt(sigma2[j]))
        }
      }
      
      # -----------------------------
      # 4. Topic regression params
      # -----------------------------
      
      delta <- rep(1.25, K)
      delta[K] <- 0
      
      Gamma <- matrix(sample(c(-1, 1), L*K, TRUE), L, K)
      Gamma[, K] <- 0
      
      # -----------------------------
      # 5. Build ψ (topic-word dist)
      # -----------------------------
      
      num_free <- length(self$V_free)
      psi_free <- matrix(0, num_free, self$K)
      
      for (k in 1:K) {
        
        start_word <- (k-1)*num_free/self$K + 1
        end_word <- k*num_free/self$K
        
        topic_words <- sample(start_word:end_word, 0.1*num_free/self$K)
        
        psi_free[topic_words, k] <- 1/length(topic_words)
      }
      
      psi <- self$build_psi(psi_free)
      
      # -----------------------------
      # 6. Generate W
      # -----------------------------
      
      W <- vector("list", N)
      
      for (i in 1:N) {
        
        theta <- self$compute_theta(FF[i, ], delta, Gamma)
        
        words <- integer(M)
        
        for (m in 1:M) {
          
          # sample topic
          z <- sample(1:K, 1, prob = theta)
          
          # sample word
          words[m] <- sample(1:V, 1, prob = psi[, z])
        }
        
        W[[i]] <- words
      }
      
      # -----------------------------
      # Return everything
      # -----------------------------
      
      list(
        Y = Y,
        W = W,
        lambda = lambda,
        
        # True parameters (useful for testing!)
        true = list(
          mu = mu,
          phi = phi,
          sigma2 = sigma2,
          FF = FF,
          delta = delta,
          Gamma = Gamma,
          psi = psi,
          rho = rho
        )
      )
    }
  )
)