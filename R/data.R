Rcpp::sourceCpp("src/utils.cpp")
source("R/tacfaSMutils.R")
source("R/tacfaSMsamplers.R")
source("R/tacfaSM.R")

N <- 500
J <- 10
K <- 5
L <- 2
V <- 30
M_min <- 60
M_max <- 90

mu <- rep(0, J)
sigma2 <- rep(1, J)
rho <- 0.3

phi <- rep(1, J-L)
phi_factor <- rep(1:2, each = J/L-1)
phi_fixed <- c(1, 6)

phi_args <- list(phi = phi,
                 phi_factor = phi_factor,
                 phi_fixed = phi_fixed)

gamma <- c(-2, -2)
gamma_factor <- c(1, 2)
gamma_fixed <- c(1, 2)
gamma_fixed_value <- 2
topic_baseline <- K

gamma_args <- list(gamma = gamma,
                   gamma_factor = gamma_factor,
                   gamma_fixed = gamma_fixed,
                   gamma_fixed_value = gamma_fixed_value,
                   topic_baseline = topic_baseline)

anchor_words <- 1:K
anchor_prob <- 0.2
nonanchor_prob <- 0.0001
num_topic_words <- 3

Psi_args <- list(anchor_words = anchor_words,
                 anchor_prob = anchor_prob,
                 nonanchor_prob = nonanchor_prob,
                 num_topic_words = num_topic_words,
                 V = V)

beta <- rep(0.01, V-K)

tau <- c(0.5, 0.5, 0.5, 0.5)
tau_full <- c(tau, 0)
tau_args <- list(tau = tau,
                 topic_baseline = topic_baseline)

data <- make_tacfaSM_data(N,
                          V,
                          M_min,
                          M_max,
                          mu,
                          sigma2,
                          rho,
                          phi_args,
                          Psi_args,
                          gamma_args,
                          tau_args)

Y <- data$Y
W <- data$W

mu_0 = 0;
sigma_mu2 = 1;
mu_phi = 0;
sigma_phi2 = 1;
a_sigma2 = 2;
b_sigma2 = 2;
a_rho = 1;
b_rho = 1;
mu_gamma = 0;
sigma_gamma2 = 1;
mu_tau = 0;
sigma_tau2 = 1;
beta_const = 0.01;
delta = 1e-05;
target_accept = 0.574;
kappa = 0.6;
niter = 30000;
nburn = 10000

fit <- tacfaSM(Y, W, V, phi_factor, phi_fixed, gamma_factor, gamma_fixed,
               gamma_fixed_value, topic_baseline, anchor_words, anchor_prob,
               nonanchor_prob, niter = 30000, nburn = 10000, warmup = 250)

