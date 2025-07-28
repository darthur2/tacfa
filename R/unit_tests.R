library(Rcpp)
library(RcppEigen)
library(Matrix)
library(MASS)
source("R/tacfa_utils.R")
source("R/mcmc_utils.R")
sourceCpp("src/tacfa_utils.cpp")

### Create Data ###

set.seed(100)

D <- 500
K <- 5
M <- sample(20:50, D, TRUE)
Z <- vector("list", D)
L <- 2
J <- 10
rho <- 0.3
a_sigma2 <- 2
b_sigma2 <- 2

R <- diag(1-rho, L) + rho
R_inv <- solve(R)
FF <- mvrnorm(D, rep(0, L), R)
Gamma <- matrix(rnorm(L*(K-1), 0, 0.5), L, K-1)
Gamma <- cbind(Gamma, 0)
Phi <- rbind(c(1, 0, abs(rnorm(J/2 - 1)), rep(0, J/2 - 1)),
             c(0, 1, rep(0, J/2 - 1), rnorm(J/2 - 1)))
delta <- rnorm(K-1, 0, 0.5)
delta <- c(delta, 0)
mu <- rnorm(J)
sigma2 <- 1/rgamma(J, a_sigma2, b_sigma2)

Theta <- exp(t(t(FF%*%Gamma) + delta))/(rowSums(exp(t(t(FF%*%Gamma) + delta))))

for (d in 1:D){
  Z[[d]] <- Matrix(t(rmultinom(M[d], 1, Theta[d,]))) 
}

Mu <- t(t(FF%*%Phi) + mu)
Sigma2 <- sapply(1:J, function(j) rnorm(D, 0, sigma2[j]))

Y <- Mu + Sigma2

### Test Z_sum Function ###

Z_sum_R <- sum_Z_R(Z)
Z_sum_cpp <- sum_Z_cpp(Z)
all.equal(Z_sum_R, Z_sum_cpp)

Z_sum <- Z_sum_cpp

system.time(replicate(100, sum_Z_R(Z)))
system.time(replicate(100, sum_Z_cpp(Z)))

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

gamma_seq <- seq(-6, 6, length = 1000)

plot(gamma_seq, sapply(gamma_seq, function(gamma_lk){
  comp_lgamma_lk_cpp(gamma_lk, FF, Gamma, Z_sum, delta, 2, 1)
  }), type = 'l')

### Test sample_gamma_lk Function ###

gamma_sd <- 0.01
ndraws <- 1000
l <- 2
k <- 4
Gamma_temp <- Gamma

gamma_lk_draws <- rep(0, ndraws)

for (iter in 1:ndraws){
  Gamma_temp[l,k] <- sample_gamma_lk(FF, Gamma_temp, Z_sum, delta, gamma_sd, l, k)
  gamma_lk_draws[iter] <- Gamma_temp[l,k]
}

### Test sample_Gamma Function ###

gamma_sd <- 0.01
ndraws <- 1000
Gamma_temp <- Gamma

Gamma_draws <- array(0, dim = c(L, K, ndraws))

for (iter in 1:ndraws){
  Gamma_temp <- sample_Gamma(FF, Gamma_temp, Z_sum, delta, gamma_sd)
  Gamma_draws[,,iter] <- Gamma_temp
}

### Test comp_lf_d Function ###

f_d_seq <- seq(-3, 3, length = 1000)
d <- 30

f_d_out <- sapply(f_d_seq, function(f_d){
  f_d <- rbind(c(f_d, FF[d,2]))
  comp_lf_d_R(f_d, Y[d,], Z_sum[d,], Phi, Gamma, R_inv, mu, sigma2, delta)
})

### Test sample_f_d Function ###

### Test sample_FF Function ##

ndraws <- 1000
sigma1 <- 0.25
sigma2 <- 0.25
f_Sigma <- rbind(c(sigma1, sigma1*sigma2*rho),
                 c(sigma1*sigma2*rho, sigma2))
FF_temp <- FF

FF_draws <- array(0, dim = c(D, L, ndraws))

for (iter in 1:ndraws){
  FF_temp <- sample_FF(FF_temp, Y, Z_sum, Phi, Gamma, R_inv, 
                       mu, sigma2, delta, f_Sigma)
  FF_draws[,,iter] <- FF_temp
}
