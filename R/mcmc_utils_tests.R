library(Rcpp)
library(RcppEigen)
library(Matrix)
library(MASS)
source("R/tacfa_utils.R")
source("R/mcmc_utils.R")
sourceCpp("src/tacfa_utils.cpp")
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
beta <- rep(0.1, V)

R <- diag(1-rho, L) + rho
R_inv <- solve(R)
FF <- mvrnorm(D, rep(0, L), R)
Gamma <- matrix(rnorm(L*(K-1), 0, 0.5), L, K-1)
Gamma <- cbind(Gamma, 0)
Phi <- rbind(c(1, 0, abs(rnorm(J/2 - 1)), rep(0, J/2 - 1)),
             c(0, 1, rep(0, J/2 - 1), rnorm(J/2 - 1)))
Psi <- rdirichlet_cpp(K, beta)
delta <- rnorm(K-1, 0, 0.5)
delta <- c(delta, 0)
mu <- rnorm(J)
sigma2 <- 1/rgamma(J, a_sigma2, b_sigma2)

Theta <- exp(t(t(FF%*%Gamma) + delta))/(rowSums(exp(t(t(FF%*%Gamma) + delta))))

for (d in 1:D){
  Z[[d]] <- sample(1:K, M[d], TRUE, Theta[d,])
  W[[d]] <- apply(Psi[,Z[[d]]], 2, function(p){
    sample(1:V, 1, FALSE, p)
  })
}

Mu <- t(t(FF%*%Phi) + mu)
Sigma2 <- sapply(1:J, function(j) rnorm(D, 0, sigma2[j]))

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
mean((Gamma_hat_R - Gamma)^2)

Gamma_hat_cpp <- apply(Gamma_draws_cpp[,,nburn:ndraws], c(1, 2), mean)
mean((Gamma_hat_cpp - Gamma)^2)

system.time(replicate(100, sample_Gamma_R(FF, Gamma_cur_R, Z_sum, delta, gamma_sd)))
system.time(replicate(100, sample_Gamma_cpp(FF, Gamma_cur_R, Z_sum, delta, gamma_sd)))

### Test sample_FF Function ###

ndraws <- 1000
nburn <- 200
f_sd <- 0.5
f_rho <- 0.5
f_Sigma <- matrix(f_sd^2*f_rho, L, L)
diag(f_Sigma) <- f_sd^2

FF_draws <- array(0, dim = c(D, L, ndraws))
FF_cur <- mvrnorm(D, rep(0, L), R)

for (iter in 1:ndraws){
  FF_cur <- sample_FF(FF_cur, Y, Z_sum, Phi, Gamma, R_inv,
                      mu, sigma2, delta, f_Sigma)
  FF_draws[,,iter] <- FF_cur
}

FF_hat <- apply(FF_draws[,,nburn:ndraws], c(1, 2), mean)
plot(FF_hat[,1], FF[,1])
plot(FF_hat[,2], FF[,2])

### Test sample_delta Function ###

ndraws <- 1000
nburn <- 200
delta_sd <- 0.1

delta_draws <- matrix(0, ndraws, K)
delta_cur <- rnorm(K-1)
delta_cur <- c(delta_cur, 0)

for (iter in 1:ndraws){
  delta_cur <- sample_delta(FF, Gamma, Z_sum, delta_cur, delta_sd)
  delta_draws[iter,] <- delta_cur
}

delta_hat <- colMeans(delta_draws[nburn:ndraws,])

mean((delta_hat - delta)^2)


