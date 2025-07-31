sum_Z_R <- function(Z, K) {
  D <- length(Z)
  Z_sum <- matrix(0, nrow = D, ncol = K)
  
  for (d in seq_len(D)) {
    z_d <- Z[[d]]
    for (k in z_d) {
      Z_sum[d,k] <- Z_sum[d,k] + 1
    }
  }
  
  Z_sum
}

rdirichlet_R <- function(n, shape){
  m <- length(shape)
  draws <- matrix(0, m, n)
  
  for (i in 1:n){
    gamma_draws_i <- rgamma(m, shape)
    draws[,i] <- gamma_draws_i/sum(gamma_draws_i)
  }
  
  draws
}

comp_WZ_counts_R <- function(W, Z, V, K) {
  WZ_counts <- matrix(0, nrow = V, ncol = K)
  
  for (d in seq_along(W)) {
    w_d <- W[[d]]
    z_d <- Z[[d]]
    for (i in seq_along(w_d)) {
      word <- w_d[i]
      topic <- z_d[i]
      WZ_counts[word, topic] <- WZ_counts[word, topic] + 1
    }
  }
  
  WZ_counts  
}
