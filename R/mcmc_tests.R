library(Rcpp)
library(RcppEigen)
library(Matrix)
source("R/mcmc_utils.R")
sourceCpp("src/mcmc_utils.cpp")

### Create Data ###

set.seed(100)

D <- 500
K <- 5
M <- sample(20:50, D, TRUE)
Z <- vector("list", D)
W <- vector("list", D)
L <- 2
J <- 10
V <- 1000
rho <- 0.3
a_sigma2 <- 2
b_sigma2 <- 2
sigma2_mu <- 4
beta <- matrix(0.1, V, K)

R <- diag(1-rho, L) + rho
R_inv <- solve(R)
FF <- rmvnorm_cpp(D, rep(0, L), R)
Gamma <- matrix(rnorm(L*(K-1), 0, 0.5), L, K-1)
Gamma <- cbind(Gamma, 0)
Phi <- rbind(c(1, 0, abs(rnorm(J/2 - 1)), rep(0, J/2 - 1)),
             c(0, 1, rep(0, J/2 - 1), abs(rnorm(J/2 - 1))))
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

Z_sum <- sum_Z_cpp(Z, K)

### Test sample_Gamma Function ###

ndraws <- 1000
nburn <- 200
gamma_sd <- 1/sqrt(D)

Gamma_draws_R <- array(0, dim = c(L, K, ndraws))
Gamma_draws_cpp <- array(0, dim = c(L, K, ndraws))
Gamma_cur_R <- matrix(rnorm(L*(K-1), 0, 0.5), L, K-1)
Gamma_cur_R <- cbind(Gamma_cur_R, 0)
Gamma_cur_cpp <- Gamma_cur_R

for (iter in 1:ndraws){
  Gamma_cur_R <- sample_Gamma_R(FF, Gamma_cur_R, Z_sum, delta, gamma_sd)
  Gamma_cur_cpp <- sample_Gamma_cpp(FF, Gamma_cur_cpp, Z_sum, delta, gamma_sd)
  Gamma_draws_R[,,iter] <- Gamma_cur_R
  Gamma_draws_cpp[,,iter] <- Gamma_cur_cpp
}

Gamma_hat_R <- apply(Gamma_draws_R[,,nburn:ndraws], c(1, 2), mean)
cor(c(Gamma_hat_R), c(Gamma))

Gamma_hat_cpp <- apply(Gamma_draws_cpp[,,nburn:ndraws], c(1, 2), mean)
cor(c(Gamma_hat_cpp), c(Gamma))

system.time(replicate(100, sample_Gamma_R(FF, Gamma_cur_R, Z_sum, delta, gamma_sd)))
system.time(replicate(100, sample_Gamma_cpp(FF, Gamma_cur_R, Z_sum, delta, gamma_sd)))

### Test sample_FF Function ###

ndraws <- 1000
nburn <- 200
f_sd <- 0.5
f_rho <- 0.5
f_Sigma <- matrix(f_sd^2*f_rho, L, L)
diag(f_Sigma) <- f_sd^2

FF_draws_R <- array(0, dim = c(D, L, ndraws))
FF_draws_cpp <- array(0, dim = c(D, L, ndraws))
FF_cur_R <- rmvnorm_cpp(D, rep(0, L), R)
FF_cur_cpp <- rmvnorm_cpp(D, rep(0, L), R)

for (iter in 1:ndraws){
  FF_cur_R <- sample_FF_R(FF_cur_R, Y, Z_sum, Phi, Gamma, R_inv,
                          mu, sigma2, delta, f_Sigma)
  FF_cur_cpp <- sample_FF_cpp(FF_cur_cpp, Y, Z_sum, Phi, Gamma, R_inv,
                              mu, sigma2, delta, f_Sigma)
  FF_draws_R[,,iter] <- FF_cur_R
  FF_draws_cpp[,,iter] <- FF_cur_cpp
}

FF_hat_R <- apply(FF_draws_R[,,nburn:ndraws], c(1, 2), mean)
FF_hat_cpp <- apply(FF_draws_cpp[,,nburn:ndraws], c(1, 2), mean)

cor(FF_hat_R[,1], FF[,1])
cor(FF_hat_R[,2], FF[,2])
cor(FF_hat_cpp[,1], FF[,1])
cor(FF_hat_cpp[,2], FF[,2])

system.time(replicate(100, sample_FF_R(FF_cur_R, Y, Z_sum, Phi, Gamma, R_inv,
                                       mu, sigma2, delta, f_Sigma)))
system.time(replicate(100, sample_FF_cpp(FF_cur_cpp, Y, Z_sum, Phi, Gamma, R_inv,
                                         mu, sigma2, delta, f_Sigma)))

### Test sample_delta Function ###

ndraws <- 1000
nburn <- 200
delta_sd <- 0.1

delta_draws_R <- matrix(0, ndraws, K)
delta_draws_cpp <- matrix(0, ndraws, K)

delta_cur_R <- rnorm(K-1)
delta_cur_R <- c(delta_cur_R, 0)
delta_cur_cpp <- rnorm(K-1)
delta_cur_cpp <- c(delta_cur_cpp, 0)

for (iter in 1:ndraws){
  delta_cur_R <- sample_delta_R(FF, Gamma, Z_sum, delta_cur_R, delta_sd)
  delta_cur_cpp <- sample_delta_cpp(FF, Gamma, Z_sum, delta_cur_cpp, delta_sd)
  delta_draws_R[iter,] <- delta_cur_R
  delta_draws_cpp[iter,] <- delta_cur_cpp
}

delta_hat_R <- colMeans(delta_draws_R[nburn:ndraws,])
delta_hat_cpp <- colMeans(delta_draws_cpp[nburn:ndraws,])

cor(delta_hat_R, delta)
cor(delta_hat_cpp, delta)

### Test sample_Z Function ###

ndraws <- 50

for (iter in 1:ndraws){
  Z_R <- sample_Z_R(W, FF, Gamma, Psi, delta)
  Z_cpp <- sample_Z_cpp(W, FF, Gamma, Psi, delta)
}

median(sapply(1:D, function(d){
  mean(Z_R[[d]] == Z[[d]])
}))

median(sapply(1:D, function(d){
  mean(Z_cpp[[d]] == Z[[d]])
}))

system.time(replicate(50, sample_Z_R(W, FF, Gamma, Psi, delta)))
system.time(replicate(50, sample_Z_cpp(W, FF, Gamma, Psi, delta)))

### Test sample_Phi Function ###

ndraws <- 1000
Phi_draws_R <- array(0, dim = c(L, J, ndraws))
Phi_draws_cpp <- array(0, dim = c(L, J, ndraws))

for (iter in 1:ndraws){
  Phi_draws_R[,,iter] <- sample_Phi_R(Y, FF, mu, sigma2, Phi, Phi_indices)
  Phi_draws_cpp[,,iter] <- sample_Phi_cpp(Y, FF, mu, sigma2, Phi, Phi_indices)
}

Phi_hat_R <- apply(Phi_draws_R, c(1, 2), mean)
Phi_hat_cpp <- apply(Phi_draws_cpp, c(1, 2), mean)

cor(c(Phi_hat_R), c(Phi))
cor(c(Phi_hat_cpp), c(Phi))

system.time(replicate(1000, sample_Phi_R(Y, FF, mu, sigma2, Phi, Phi_indices)))
system.time(replicate(1000, sample_Phi_cpp(Y, FF, mu, sigma2, Phi, Phi_indices)))

### Test sample_sigma2 Function ###

ndraws <- 1000
sigma2_draws_R <- matrix(0, ndraws, J)
sigma2_draws_cpp <- matrix(0, ndraws, J)

for (iter in 1:ndraws){
  sigma2_draws_R[iter,] <- sample_sigma2_R(Y, mu, FF, Phi, a_sigma2, b_sigma2)
  sigma2_draws_cpp[iter,] <- sample_sigma2_cpp(Y, mu, FF, Phi, a_sigma2, b_sigma2)
}

sigma2_hat_R <- colMeans(sigma2_draws_R)
sigma2_hat_cpp <- colMeans(sigma2_draws_cpp)

cor(sigma2_hat_R, sigma2)
cor(sigma2_hat_cpp, sigma2)

system.time(replicate(1000, sample_sigma2_R(Y, mu, FF, Phi, a_sigma2, b_sigma2)))
system.time(replicate(1000, sample_sigma2_cpp(Y, mu, FF, Phi, a_sigma2, b_sigma2)))

### Test sample_mu Function ###

ndraws <- 1000
mu_draws_R <- matrix(0, ndraws, J)
mu_draws_cpp <- matrix(0, ndraws, J)

for (iter in 1:ndraws){
  mu_draws_R[iter,] <- sample_mu_R(Y, sigma2, FF, Phi, sigma2_mu)
  mu_draws_cpp[iter,] <- sample_mu_cpp(Y, sigma2, FF, Phi, sigma2_mu)
}

mu_hat_R <- colMeans(mu_draws_R)
mu_hat_cpp <- colMeans(mu_draws_cpp)

cor(mu_hat_R, mu)
cor(mu_hat_cpp, mu)

system.time(replicate(1000, sample_mu_R(Y, sigma2, FF, Phi, sigma2_mu)))
system.time(replicate(1000, sample_mu_cpp(Y, sigma2, FF, Phi, sigma2_mu)))

### Test sample_Psi Function ###

ndraws <- 1000
Psi_draws_R <- array(0, dim = c(V, K, ndraws))
Psi_draws_cpp <- array(0, dim = c(V, K, ndraws))

for (iter in 1:ndraws){
  Psi_draws_R[,,iter] <- sample_Psi_R(W, Z, beta)
  Psi_draws_cpp[,,iter] <- sample_Psi_cpp(W, Z, beta)
}

Psi_hat_R <- apply(Psi_draws_R, c(1, 2), mean)
Psi_hat_cpp <- apply(Psi_draws_cpp, c(1, 2), mean)

diag(cor(Psi_hat_R, Psi))
diag(cor(Psi_hat_cpp, Psi))

system.time(replicate(1000, sample_Psi_R(W, Z, beta)))
system.time(replicate(1000, sample_Psi_cpp(W, Z, beta)))
