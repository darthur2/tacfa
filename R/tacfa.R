library(Rcpp)
library(RcppEigen)
library(Matrix)
source("R/mcmc_utils.R")
sourceCpp("src/mcmc_utils.cpp")

tacfa <- function(Y, W, Z, V, L, K, Phi_indices,
                  init = NULL, prior = NULL, mh = NULL,
                  niter = 2000, nburn = 200){
  D <- nrow(Y)
  J <- ncol(Y)
  
  FF_draws <- array(0, dim = c(D, L, niter-nburn))
  Gamma_draws <- array(0, dim = c(L, K, niter-nburn))
  delta_draws <- matrix(0, niter-nburn, K)
  Psi_draws <- array(0, dim = c(V, K, niter-nburn))
  Phi_draws <- array(0, dim = c(L, J, niter-nburn))
  mu_draws <- matrix(0, niter-nburn, J)
  sigma2_draws <- matrix(0, niter-nburn, J)
  
  if (is.null(prior)){
    beta <- matrix(0.1, V, K)
    a_sigma2 <- b_sigma2 <- 2
    rho <- 0.5
    R <- diag(1-rho, L) + rho
    R_inv <- solve(R)
  }
  
  if (is.null(init)){
    Gamma <- cbind(diag(L), matrix(rnorm(L*(K-L)), L, K-L))
    delta <- rnorm(K-1)
    delta <- c(delta, 0)
    Psi <- rdirichlet_cpp(K, beta)
    Phi <- rbind(c(1, 0, rep(c(1, 0), each = (J-L)/2)),
                 c(0, 1, rep(c(0, 1), each = (J-L)/2)))
    mu <- colMeans(Y)
    sigma2 <- apply(Y, 2, var)
    FF <- sapply(1:L, function(l) colSums((t(Y)-mu))/J)
  }
  
  if (is.null(mh)){
    f_sd <- 1
    f_rho <- 0.5
    f_Sigma <- matrix(f_sd^2*f_rho, L, L)
    diag(f_Sigma) <- f_sd^2
    # f_Sigma <- 1/(J+K-1)*f_Sigma
    f_Sigma <- 1/sqrt(D)*f_Sigma
    
    gamma_sd <- sqrt(1/(D + 1))
    delta_sd <- sqrt(1/(D + 1))
  }
  
  # Z <- vector("list", D)
  # Theta <- exp(FF%*%Gamma)/rowSums(exp(FF%*%Gamma))
  # 
  # for (d in 1:D){
  #   Z[[d]] <- sample(1:K, length(W[[d]]), TRUE, Theta[d,])
  # }
  
  for (iter in 1:niter){
    Z_sum <- sum_Z_cpp(Z, K)
    Gamma <- sample_Gamma_cpp(FF, Gamma, Z_sum, delta, gamma_sd)
    Gamma[1:L, 1:L] <- diag(L)
    Gamma[2,3:4] <- 0
    Gamma[1,5:6] <- 0
    delta <- sample_delta_cpp(FF, Gamma, Z_sum, delta, delta_sd)
    Psi <- sample_Psi_cpp(W, Z, beta)
    # Z <- sample_Z_cpp(W, FF, Gamma, Psi, delta)
    FF <- sample_FF_cpp(FF, Y, Z_sum, Phi, Gamma, R_inv, mu, sigma2, delta, f_Sigma)
    Phi <- sample_Phi_cpp(Y, FF, mu, sigma2, Phi, Phi_indices)
    mu <- sample_mu_cpp(Y, sigma2, FF, Phi, sigma2_mu)
    sigma2 <- sample_sigma2_cpp(Y, mu, FF, Phi, a_sigma2, b_sigma2)
    
    if (iter > nburn){
      FF_draws[,,iter-nburn] <- FF
      Gamma_draws[,,iter-nburn] <- Gamma
      delta_draws[iter-nburn,] <- delta
      Psi_draws[,,iter-nburn] <- Psi
      Phi_draws[,,iter-nburn] <- Phi
      mu_draws[iter-nburn,] <- mu
      sigma2_draws[iter-nburn,] <- sigma2
    }
  }
  
  list(FF_draws = FF_draws,
       Gamma_draws = Gamma_draws,
       delta_draws = delta_draws,
       Psi_draws = Psi_draws,
       Phi_draws = Phi_draws,
       mu_draws = mu_draws,
       sigma2_draws = sigma2_draws)
}

### Create Data ###

D <- 750
K <- 6
M <- sample(150:200, D, TRUE)
Z <- vector("list", D)
W <- vector("list", D)
L <- 2
J <- 10
V <- 100
rho <- 0.3
a_sigma2 <- 2
b_sigma2 <- 2
sigma2_mu <- 4
beta <- matrix(0.1, V, K)

R <- diag(1-rho, L) + rho
R_inv <- solve(R)
FF <- rmvnorm_cpp(D, rep(0, L), R)
# Gamma <- cbind(diag(L), matrix(rnorm(L*(K-L)), L, K-L))
Gamma <- rbind(c(1, 0, rep(c(1, 0), each = floor((K-L)/2))),
               c(0, 1, rep(c(0, 1), each = ceiling((K-L)/2))))
# Phi <- rbind(c(1, 0, abs(rnorm(J/2 - 1)), rep(0, J/2 - 1)),
#              c(0, 1, rep(0, J/2 - 1), abs(rnorm(J/2 - 1))))
Phi <- rbind(c(1, 0, rep(c(1, 0), each = (J-L)/2)),
             c(0, 1, rep(c(0, 1), each = (J-L)/2)))
Phi_indices <- list(3:6, 7:10) # Generalize later
Psi <- matrix(0, V, K)

for (k in 1:K){
  Psi[,k] <- rdirichlet_cpp(1, beta[,k])  
}

delta <- rnorm(K-1, 0, 0.5)
delta <- c(delta, 0)

mu <- rnorm(J, 0, sqrt(sigma2_mu))
sigma2 <- 1/rgamma(J, a_sigma2, b_sigma2)

Theta <- exp(t(t(FF%*%Gamma) + delta))/(rowSums(exp(t(t(FF%*%Gamma) + delta))))

for (d in 1:D){
  Z[[d]] <- sample(1:K, M[d], TRUE, Theta[d,])
  W[[d]] <- apply(Psi[,Z[[d]]], 2, function(p){
    sample(1:V, 1, FALSE, p)
  })
}

Mu <- t(t(FF%*%Phi) + mu)
Sigma2 <- sapply(1:J, function(j) rnorm(D, 0, sqrt(sigma2[j])))

Y <- Mu + Sigma2

### Fit Model ###

tacfa_fit <- tacfa(Y, W, Z, V, L, K, Phi_indices, niter = 2000, nburn = 200)

### Analyze Results ###

FF_draws <- tacfa_fit$FF_draws
Gamma_draws <- tacfa_fit$Gamma_draws
delta_draws <- tacfa_fit$delta_draws
Psi_draws <- tacfa_fit$Psi_draws
Phi_draws <- tacfa_fit$Phi_draws
mu_draws <- tacfa_fit$mu_draws
sigma2_draws <- tacfa_fit$sigma2_draws

D <- nrow(Y)
J <- ncol(Y)

FF_draws <- array(0, dim = c(D, L, niter-nburn))
Gamma_draws <- array(0, dim = c(L, K, niter-nburn))
delta_draws <- matrix(0, niter-nburn, K)
Psi_draws <- array(0, dim = c(V, K, niter-nburn))
Phi_draws <- array(0, dim = c(L, J, niter-nburn))
mu_draws <- matrix(0, niter-nburn, J)
sigma2_draws <- matrix(0, niter-nburn, J)

if (is.null(prior)){
  beta <- matrix(0.1, V, K)
  a_sigma2 <- b_sigma2 <- 2
  rho <- 0.5
  R <- diag(1-rho, L) + rho
  R_inv <- solve(R)
}

if (is.null(init)){
  Gamma <- cbind(diag(L), matrix(rnorm(L*(K-L)), L, K-L))
  delta <- rnorm(K-1)
  delta <- c(delta, 0)
  Psi <- rdirichlet_cpp(K, beta)
  Phi <- rbind(c(1, 0, rep(c(1, 0), each = (J-L)/2)),
               c(0, 1, rep(c(0, 1), each = (J-L)/2)))
  mu <- colMeans(Y)
  sigma2 <- apply(Y, 2, var)
  FF <- sapply(1:L, function(l) colSums((t(Y)-mu))/J)
}

if (is.null(mh)){
  f_sd <- 1
  f_rho <- 0.5
  f_Sigma <- matrix(f_sd^2*f_rho, L, L)
  diag(f_Sigma) <- f_sd^2
  # f_Sigma <- 1/(J+K-1)*f_Sigma
  f_Sigma <- 1/sqrt(D)*f_Sigma
  
  gamma_sd <- sqrt(1/(D + 1))
  delta_sd <- sqrt(1/(D + 1))
}

# Z <- vector("list", D)
# Theta <- exp(FF%*%Gamma)/rowSums(exp(FF%*%Gamma))
# 
# for (d in 1:D){
#   Z[[d]] <- sample(1:K, length(W[[d]]), TRUE, Theta[d,])
# }

for (iter in 1:niter){
  Z_sum <- sum_Z_cpp(Z, K)
  Gamma <- sample_Gamma_cpp(FF, Gamma, Z_sum, delta, gamma_sd)
  Gamma[1:L, 1:L] <- diag(L)
  Gamma[2,3:4] <- 0
  Gamma[1,5:6] <- 0
  delta <- sample_delta_cpp(FF, Gamma, Z_sum, delta, delta_sd)
  Psi <- sample_Psi_cpp(W, Z, beta)
  # Z <- sample_Z_cpp(W, FF, Gamma, Psi, delta)
  FF <- sample_FF_cpp(FF, Y, Z_sum, Phi, Gamma, R_inv, mu, sigma2, delta, f_Sigma)
  Phi <- sample_Phi_cpp(Y, FF, mu, sigma2, Phi, Phi_indices)
  
  if (any(is.na(Phi)) | any(is.infinite(Phi))){
    break
  }
  
  mu <- sample_mu_cpp(Y, sigma2, FF, Phi, sigma2_mu)
  
  if (any(is.na(mu))){
    break
  }
  
  sigma2 <- sample_sigma2_cpp(Y, mu, FF, Phi, a_sigma2, b_sigma2)
  
  if (any(is.na(mu))){
    break
  }
  
  if (iter > nburn){
    FF_draws[,,iter-nburn] <- FF
    Gamma_draws[,,iter-nburn] <- Gamma
    delta_draws[iter-nburn,] <- delta
    Psi_draws[,,iter-nburn] <- Psi
    Phi_draws[,,iter-nburn] <- Phi
    mu_draws[iter-nburn,] <- mu
    sigma2_draws[iter-nburn,] <- sigma2
  }
}

FF_acc_rates <- sapply(1:D, function(d) mean(diff(FF_draws[d,1,]) != 0))
Gamma_acc_rates <- matrix(0, L, K)

for (l in 1:L){
  for (k in 1:K){
    Gamma_acc_rates[l,k] <- mean(diff(Gamma_draws[l,k,]) != 0)
  }
}

delta_acc_rates <- sapply(1:K, function(k) mean(diff(delta_draws[,k]) != 0))

hist(FF_acc_rates)
Gamma_acc_rates
delta_acc_rates

FF_hat <- apply(FF_draws, c(1, 2), mean, na.rm = TRUE)
Gamma_hat <- apply(Gamma_draws, c(1, 2), mean, na.rm = TRUE)
delta_hat <- colMeans(delta_draws, na.rm = TRUE)
Psi_hat <- apply(Psi_draws, c(1, 2), mean, na.rm = TRUE)
Phi_hat <- apply(Phi_draws, c(1, 2), mean, na.rm = TRUE)
mu_hat <- colMeans(mu_draws, na.rm = TRUE)
sigma2_hat <- colMeans(sigma2_draws, na.rm = TRUE)

cor(FF_hat[,1], FF[,1])
cor(FF_hat[,2], FF[,2])
cor(c(Gamma_hat),c(Gamma))
cor(delta_hat, delta)
cor(Psi_hat, Psi)
cor(mu_hat, mu)
cor(sigma2_hat, sigma2)
cor(c(Phi_hat), c(Phi))
