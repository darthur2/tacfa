TacfaSMEngine <- R6Class(
  "TacfaSMEngine",
  
  public = list(
    
    model = NULL,
    state = NULL,
    
    setup = function(model, Y, W, lambda) {
      self$model <- model
      self$state <- private$initialize_state(model, Y, W, lambda)
    },
    
    run = function() {
      stop("Not implemented")
    },
    
    get_state = function() {
      self$state
    }
  ),
  
  private = list(
    
    initialize_state = function(model, Y, W, lambda) {
      
      N <- nrow(Y)
      J <- ncol(Y)
      
      state <- list(
        Y = Y,
        W = W,
        lambda = lambda,
        
        mu = rep(0, J),
        sigma2 = rep(1, J),
        phi = rep(1, J),
        
        FF = matrix(rnorm(N * model$L), N, model$L),
        
        delta = rep(0, model$K),
        Gamma = matrix(rnorm(model$L * model$K), model$L, model$K),
        
        psi_free = matrix(
          1 / length(model$V_free),
          length(model$V_free),
          model$K
        ),
        
        rho = 0.5
      )
      
      state
    }
  )
)

TacfaSMGibbs <- R6::R6Class(
  "TacfaSMGibbs",
  
  inherit = TacfaSMEngine,
  
  public = list(
    
    n_iter = NULL,
    burn_in = NULL,
    thin = NULL,
    verbose = NULL,

    samples = NULL,
    iter = NULL,
    
    adapt = NULL,
    
    initialize = function(
      n_iter = 1000,
      burn_in = 200,
      thin = 5,
      verbose = TRUE
    ) {
      self$n_iter <- n_iter
      self$burn_in <- burn_in
      self$thin <- thin
      self$verbose <- verbose
    },
    
    setup = function(model, Y, W, lambda) {
      super$setup(model, Y, W, lambda)
      
      # Initialize adaptive proposal scales
      self$adapt <- list(
        F_sd = rep(0.1, model$L),
        delta_sd = rep(0.1, model$K),
        gamma_sd = rep(0.1, model$L),
        rho_sd = 0.05
      )
    },
    
    run = function() {
      
      n_save <- floor((self$n_iter - self$burn_in) / self$thin)
      
      self$samples <- NULL
      save_idx <- 1
      
      for (iter in 1:self$n_iter) {
        
        self$iter <- iter
        
        # -----------------
        # Updates
        # -----------------
        self$update_mu()
        self$update_sigma2()
        self$update_phi()
        self$update_psi()
        
        self$update_F()
        self$update_delta()
        self$update_Gamma()
        self$update_rho()
        
        # -----------------
        # Store samples
        # -----------------
        if (iter > self$burn_in && (iter - self$burn_in) %% self$thin == 0) {
          
          flat <- private$flatten_state(self$state)
          
          # Initialize matrix on first save
          if (is.null(self$samples)) {
            self$samples <- matrix(NA, nrow = n_save, ncol = length(flat))
            colnames(self$samples) <- names(flat)
          }
          
          self$samples[save_idx, ] <- flat
          save_idx <- save_idx + 1
        }
        
        if (self$verbose && iter %% 100 == 0) {
        cat("Iteration:", iter, "\n")
        }
      }
    },
    
    get_samples = function() {
      self$samples
    },
    
    # =====================================================
    # 🔁 Gibbs updates
    # =====================================================
    
    update_mu = function() {
      
      Y <- self$state$Y
      FF <- self$state$FF
      phi <- self$state$phi
      sigma2 <- self$state$sigma2
      lambda <- self$state$lambda
      
      mu0 <- self$model$mu0
      sigma_mu2 <- self$model$sigma_mu2
      
      N <- nrow(Y)
      J <- ncol(Y)
      
      for (j in 1:J) {
        
        V_j <- 1 / (N / sigma2[j] + 1 / sigma_mu2)
        
        resid_sum <- sum(Y[, j] - phi[j] * FF[, lambda[j]])
        
        m_j <- V_j * (resid_sum / sigma2[j] + mu0 / sigma_mu2)
        
        self$state$mu[j] <- rnorm(1, m_j, sqrt(V_j))
      }
    },
    
    update_sigma2 = function() {
      
      Y <- self$state$Y
      FF <- self$state$FF
      mu <- self$state$mu
      phi <- self$state$phi
      lambda <- self$state$lambda
      
      a0 <- self$model$a_sigma2
      b0 <- self$model$b_sigma2
      
      N <- nrow(Y)
      J <- ncol(Y)
      
      for (j in 1:J) {
        
        mean_vec <- mu[j] + phi[j] * FF[, lambda[j]]
        ss <- sum((Y[, j] - mean_vec)^2)
        
        shape <- a0 + N / 2
        rate <- b0 + 0.5 * ss
        
        self$state$sigma2[j] <- private$rinvgamma(1, shape, rate)
      }
    },
    
    update_phi = function() {
      
      Y <- self$state$Y
      FF <- self$state$FF
      mu <- self$state$mu
      sigma2 <- self$state$sigma2
      lambda <- self$state$lambda
      
      mu_phi <- self$model$mu_phi
      sigma_phi2 <- self$model$sigma_phi2
      
      J <- ncol(Y)
      
      for (j in 1:J) {
        
        f_j <- FF[, lambda[j]]
        
        V_j <- 1 / (sum(f_j^2) / sigma2[j] + 1 / sigma_phi2)
        
        m_j <- V_j * (
          sum(f_j * (Y[, j] - mu[j])) / sigma2[j] +
          mu_phi / sigma_phi2
        )
        
        draw <- private$rtruncnorm_lower(1, m_j, sqrt(V_j))
        
        self$state$phi[j] <- draw
      }
    },
    
    update_psi = function() {
      
      W <- self$state$W
      FF <- self$state$FF
      delta <- self$state$delta
      Gamma <- self$state$Gamma
      
      model <- self$model
      
      V_free <- model$V_free
      K <- model$K
      
      counts <- matrix(0, length(V_free), K)
      psi <- model$build_psi(self$state$psi_free)
      
      for (i in seq_along(W)) {
        
        theta <- model$compute_theta(FF[i, ], delta, Gamma)
        
        for (w in W[[i]]) {
          
          probs <- psi[w, ] * theta
          probs <- probs / sum(probs)
          
          if (w %in% V_free) {
            idx <- match(w, V_free)
            counts[idx, ] <- counts[idx, ] + probs
          }
        }
      }
      
      for (k in 1:K) {
        
        alpha <- model$beta_V + counts[, k]
        draw <- rgamma(length(alpha), shape = alpha)
        self$state$psi_free[, k] <- draw / sum(draw)
      }
    },
    
    # =====================================================
    # 🔁 MH updates
    # =====================================================
    
    update_F = function() {
      
      N <- nrow(self$state$FF)
      L <- self$model$L
      
      for (i in 1:N) {
        
        current <- self$state$FF[i, ]
        proposal <- current + rnorm(L, 0, self$adapt$F_sd)
        
        log_acc <- private$log_post_f_i(i, proposal) -
                   private$log_post_f_i(i, current)
        
        accept <- 0
        if (log(runif(1)) < log_acc) {
          self$state$FF[i, ] <- proposal
          accept <- 1
        }
        
        if (self$iter <= self$burn_in) {
          self$adapt$F_sd <- self$adapt$F_sd *
            exp(0.01 * (accept - 0.234))
        }
      }
    },
    
    update_delta = function() {
      
      current <- self$state$delta
      proposal <- current + rnorm(length(current), 0, self$adapt$delta_sd)
      
      proposal[length(proposal)] <- 0
      
      log_acc <- private$log_post_delta(proposal) -
                 private$log_post_delta(current)
      
      accept <- 0
      if (log(runif(1)) < log_acc) {
        self$state$delta <- proposal
        accept <- 1
      }
      
      if (self$iter <= self$burn_in) {
        self$adapt$delta_sd <- self$adapt$delta_sd *
          exp(0.01 * (accept - 0.234))
      }
    },
    
    update_Gamma = function() {
      
      current <- self$state$Gamma
      proposal <- current + matrix(
        rnorm(length(current), 0, self$adapt$gamma_sd),
        nrow = nrow(current)
      )
      
      proposal[, ncol(proposal)] <- 0
      
      log_acc <- private$log_post_gamma(proposal) -
                 private$log_post_gamma(current)
      
      accept <- 0
      if (log(runif(1)) < log_acc) {
        self$state$Gamma <- proposal
        accept <- 1
      }
      
      if (self$iter <= self$burn_in) {
        self$adapt$gamma_sd <- self$adapt$gamma_sd *
          exp(0.01 * (accept - 0.234))
      }
    },
    
    update_rho = function() {
      
      current <- self$state$rho
      proposal <- current + rnorm(1, 0, self$adapt$rho_sd)
      
      if (proposal <= 0 || proposal >= 1) return()
      
      log_acc <- private$log_post_rho(proposal) -
                 private$log_post_rho(current)
      
      accept <- 0
      if (log(runif(1)) < log_acc) {
        self$state$rho <- proposal
        accept <- 1
      }
      
      if (self$iter <= self$burn_in) {
        self$adapt$rho_sd <- self$adapt$rho_sd *
          exp(0.01 * (accept - 0.44))
      }
    }
  ),
  
  private = list(
    
    rtruncnorm_lower = function(n, mean, sd, lower = 0) {
      alpha <- (lower - mean) / sd
      u <- runif(n, pnorm(alpha), 1)
      mean + sd * qnorm(u)
    },
    
    flatten_state = function(state) {
      
      out <- c()
      names_out <- c()
      
      add_vec <- function(x, name) {
        out <<- c(out, x)
        names_out <<- c(names_out, paste0(name, "[", seq_along(x), "]"))
      }
      
      add_mat <- function(x, name) {
        idx <- which(matrix(TRUE, nrow(x), ncol(x)), arr.ind = TRUE)
        out <<- c(out, as.vector(x))
        names_out <<- c(
          names_out,
          paste0(name, "[", idx[,1], ",", idx[,2], "]")
        )
      }
      
      # Scalars / vectors
      add_vec(state$mu, "mu")
      add_vec(state$sigma2, "sigma2")
      add_vec(state$phi, "phi")
      add_vec(state$delta, "delta")
      
      # Matrices
      add_mat(state$Gamma, "Gamma")
      add_mat(state$FF, "F")
      
      # Scalar
      out <- c(out, state$rho)
      names_out <- c(names_out, "rho")
      
      names(out) <- names_out
      
      return(out)
    },
    
    rinvgamma = function(n, shape, rate) {
      1 / rgamma(n, shape = shape, rate = rate)
    },
    
    log_post_f_i = function(i, f_i) {
      
      state <- self$state
      model <- self$model
      
      R <- matrix(state$rho, model$L, model$L)
      diag(R) <- 1
      Rinv <- solve(R)
      
      lp <- -0.5 * t(f_i) %*% Rinv %*% f_i
      
      for (j in 1:ncol(state$Y)) {
        mean_ij <- state$mu[j] +
          state$phi[j] * f_i[state$lambda[j]]
        
        lp <- lp - 0.5 / state$sigma2[j] *
          (state$Y[i, j] - mean_ij)^2
      }
      
      theta <- model$compute_theta(f_i, state$delta, state$Gamma)
      psi <- model$build_psi(state$psi_free)
      
      for (w in state$W[[i]]) {
        lp <- lp + log(sum(psi[w, ] * theta))
      }
      
      return(lp)
    },
    
    log_post_delta = function(delta) {
      
      state <- self$state
      model <- self$model
      
      lp <- 0
      
      psi <- model$build_psi(state$psi_free)
      
      for (i in 1:nrow(state$FF)) {
        theta <- model$softmax(delta + state$FF[i, ] %*% state$Gamma)
        
        for (w in state$W[[i]]) {
          lp <- lp + log(sum(psi[w, ] * theta))
        }
      }
      
      lp <- lp - sum((delta - model$mu_delta)^2) /
        (2 * model$sigma_delta2)
      
      return(lp)
    },
    
    log_post_gamma = function(Gamma) {
      
      state <- self$state
      model <- self$model
      
      lp <- 0
      
      psi <- model$build_psi(state$psi_free)
      
      for (i in 1:nrow(state$FF)) {
        theta <- model$softmax(state$delta + state$FF[i, ] %*% Gamma)
        
        for (w in state$W[[i]]) {
          lp <- lp + log(sum(psi[w, ] * theta))
        }
      }
      
      lp <- lp - sum((Gamma - model$mu_gamma)^2) /
        (2 * model$sigma_gamma2)
      
      return(lp)
    },
    
    log_post_rho = function(rho) {
      
      if (rho <= 0 || rho >= 1) return(-Inf)
      
      state <- self$state
      model <- self$model
      
      lp <- dbeta(rho, model$a_rho, model$b_rho, log = TRUE)
      
      L <- model$L
      R <- matrix(rho, L, L)
      diag(R) <- 1
      
      Rinv <- solve(R)
      logdet <- as.numeric(determinant(R, log = TRUE)$modulus)
      
      for (i in 1:nrow(state$FF)) {
        fi <- state$FF[i, ]
        lp <- lp - 0.5 * (t(fi) %*% Rinv %*% fi) - 0.5 * logdet
      }
      
      return(lp)
    }
  )
)