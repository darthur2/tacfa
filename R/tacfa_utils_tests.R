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

### Test Z_sum Function ###

Z_sum_R <- sum_Z_R(Z, K)
Z_sum_cpp <- sum_Z_cpp(Z, K)
all.equal(Z_sum_R, Z_sum_cpp)

Z_sum <- Z_sum_cpp

system.time(replicate(1000, sum_Z_R(Z, K)))
system.time(replicate(1000, sum_Z_cpp(Z, K)))

### Test comp_lgamma_lk Function ###

log_Gamma_R <- matrix(0, L, K)
log_Gamma_cpp <- matrix(0, L, K)

for (l in 1:L){
  for (k in 1:K){
    log_Gamma_R[l,k] <- comp_lgamma_lk_R(Gamma[l,k], FF, Gamma, Z_sum, delta, l, k)
    log_Gamma_cpp[l,k] <- comp_lgamma_lk_cpp(Gamma[l,k], FF, Gamma, Z_sum, delta, l, k)
  }
}

all.equal(log_Gamma_R, log_Gamma_cpp)

system.time(replicate(1000, comp_lgamma_lk_R(Gamma[1,1], FF, Gamma, Z_sum, delta, l, k)))
system.time(replicate(1000, comp_lgamma_lk_cpp(Gamma[1,1], FF, Gamma, Z_sum, delta, l, k)))

### Test comp_lf_d_R ###

log_FF_R <- matrix(0, D, L)
log_FF_cpp <- matrix(0, D, L)

for (d in 1:D){
  for (l in 1:L){
    log_FF_R[d,l] <- comp_lf_d_R(FF[d,], Y[d,], Z_sum[d,], Phi, Gamma,
                                 R_inv, mu, sigma2, delta)
    log_FF_cpp[d,l] <- comp_lf_d_cpp(FF[d,], Y[d,], Z_sum[d,], Phi, Gamma,
                                     R_inv, mu, sigma2, delta)
  }
}

all.equal(log_FF_R, log_FF_cpp)

system.time(replicate(10000, comp_lf_d_R(FF[d,], Y[d,], Z_sum[d,], Phi, Gamma,
                                         R_inv, mu, sigma2, delta)))
system.time(replicate(10000, comp_lf_d_cpp(FF[d,], Y[d,], Z_sum[d,], Phi, Gamma,
                                           R_inv, mu, sigma2, delta)))

### Test comp_ldelta_k Function ###

delta_R <- rep(0, K)
delta_cpp <- rep(0, K)

for (k in 1:K){
  delta_R[k] <- comp_ldelta_k_R(delta[k], FF, Gamma, Z_sum, delta, k)
  delta_cpp[k] <- comp_ldelta_k_cpp(delta[k], FF, Gamma, Z_sum, delta, k)
}

all.equal(delta_R, delta_cpp)

system.time(replicate(1000, comp_ldelta_k_R(delta[k], FF, Gamma, Z_sum, delta, k)))
system.time(replicate(1000, comp_ldelta_k_cpp(delta[k], FF, Gamma, Z_sum, delta, k)))

### Test rdirchlet Function ###

shape <- rep(0.1, V)
dir_draws_R <- rdirichlet_R(K, shape)
dir_draws_cpp <- rdirichlet_cpp(K, shape)

system.time(replicate(1000, rdirichlet_R(K, shape)))
system.time(replicate(1000, rdirichlet_cpp(K, shape)))

### Test comp_Z_post Function ###

Z_post_R <- comp_Z_post_R(W, FF, Gamma, Psi, delta) 
Z_post_cpp <- comp_Z_post_cpp(W, FF, Gamma, Psi, delta)

all.equal(Z_post_R, Z_post_cpp)

system.time(replicate(100, comp_Z_post_R(W, FF, Gamma, Psi, delta)))
system.time(replicate(100, comp_Z_post_cpp(W, FF, Gamma, Psi, delta)))

### Test comp_WZ_counts Function ###

WZ_counts_R <- comp_WZ_counts_R(W, Z, V, K)
WZ_counts_cpp <- comp_WZ_counts_cpp(W, Z, V, K)

all.equal(WZ_counts_R, WZ_counts_cpp)

system.time(replicate(1000, comp_WZ_counts_R(W, Z, V, K)))
system.time(replicate(1000, comp_WZ_counts_cpp(W, Z, V, K)))
