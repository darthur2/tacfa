#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::VectorXd sample_mu_cpp(
    const Eigen::MatrixXd& Y,
    const Eigen::MatrixXd& FF,
    const Eigen::MatrixXd& Phi,
    const Eigen::VectorXd& sigma2,
    double mu_0,
    double sigma_mu2
) {
  
  int N = Y.rows();
  int J = Y.cols();
  
  // Compute FF %*% Phi
  Eigen::MatrixXd FF_Phi = FF * Phi;
  
  // Column sums of residuals
  Eigen::VectorXd colsum = (Y - FF_Phi).colwise().sum();
  
  Eigen::VectorXd mu(J);
  
  for (int j = 0; j < J; ++j) {
    
    double mu_var_j = 1.0 / (N / sigma2[j] + 1.0 / sigma_mu2);
    
    double mu_mean_j = mu_var_j * (
      (1.0 / sigma2[j]) * colsum[j] +
        mu_0 / sigma_mu2
    );
    
    double sd_j = std::sqrt(mu_var_j);
    
    mu[j] = R::rnorm(mu_mean_j, sd_j);
  }
  
  return mu;
}

// [[Rcpp::export]]
Eigen::VectorXd sample_sigma2_cpp(
    const Eigen::MatrixXd& Y,
    const Eigen::MatrixXd& FF,
    const Eigen::MatrixXd& Phi,
    const Eigen::VectorXd& mu,   // length J
    double a_sigma2,
    double b_sigma2
) {
  
  int N = Y.rows();
  int J = Y.cols();
  
  // Compute FF %*% Phi once
  Eigen::MatrixXd FF_Phi = FF * Phi;
  
  // Residuals with broadcasting (no Mu matrix created)
  Eigen::MatrixXd resid = Y - FF_Phi;
  resid.rowwise() -= mu.transpose();  // subtract mu from each row
  
  // Column-wise sum of squares
  Eigen::VectorXd resid_sum2 = resid.array().square().colwise().sum();
  
  // Posterior parameters
  double sigma2_shape = a_sigma2 + N / 2.0;
  Eigen::VectorXd sigma2_rate = (0.5 * resid_sum2.array() + b_sigma2).matrix();
  
  // Sample sigma2
  Eigen::VectorXd sigma2(J);
  
  for (int j = 0; j < J; ++j) {
    double scale = 1.0 / sigma2_rate[j];  // convert rate → scale
    sigma2[j] = 1.0 / R::rgamma(sigma2_shape, scale);
  }
  
  return sigma2;
}

// helper: truncated normal (lower = 0)
double rtnorm_lower_single_cpp(double mu, double sd, double lower) {
  double x;
  do {
    x = R::rnorm(mu, sd);
  } while (x <= lower);
  return x;
}

// [[Rcpp::export]]
Eigen::VectorXd sample_phi_cpp(
    const Eigen::MatrixXd& Y,
    const Eigen::MatrixXd& FF,
    const Eigen::VectorXd& mu,
    const Eigen::VectorXd& sigma2,
    double mu_phi,
    double sigma_phi2,
    const IntegerVector& phi_fixed,   // 1-based
    const IntegerVector& phi_factor   // length J-L, 1-based
) {
  
  int N = Y.rows();
  int J = Y.cols();
  int L = FF.cols();
  
  // ---- Safety checks ----
  if (mu.size() != J)
    Rcpp::stop("mu must have length J");
  
  if (sigma2.size() != J)
    Rcpp::stop("sigma2 must have length J");
  
  // ---- Initialize phi ----
  Eigen::VectorXd phi = Eigen::VectorXd::Ones(J-L);
  
  // ---- Residuals ----
  Eigen::MatrixXd resid = Y;
  resid.rowwise() -= mu.transpose();
  
  // ---- Mark fixed indices ----
  std::vector<bool> is_fixed(J, false);
  
  for (int i = 0; i < phi_fixed.size(); ++i) {
    int idx = phi_fixed[i] - 1;
    if (idx < 0 || idx >= J) {
      Rcpp::stop("phi_fixed index out of bounds");
    }
    is_fixed[idx] = true;
  }
  
  // ---- Precompute cross-products ----
  Eigen::MatrixXd XtX = FF.transpose() * FF;     // L × L
  Eigen::MatrixXd XtR = FF.transpose() * resid;  // L × J
  
  // ---- Main loop (only scalar work now) ----
  int free_idx = 0;
  
  for (int j = 0; j < J; ++j) {
    
    if (is_fixed[j]) continue;
    
    if (free_idx >= phi_factor.size()) {
      Rcpp::stop("phi_factor too short for number of free phi");
    }
    
    int k = phi_factor[free_idx] - 1;
    
    if (k < 0 || k >= L) {
      Rcpp::stop("Invalid phi_factor index at free_idx = %d", free_idx);
    }
    
    // ---- Sufficient statistics (fast lookup) ----
    double sum_FF2 = XtX(k, k);
    double sum_FF_resid = XtR(k, j);
    
    // ---- Posterior variance ----
    double phi_var = 1.0 / (
      (1.0 / sigma2[j]) * sum_FF2 +
        1.0 / sigma_phi2
    );
    
    if (!std::isfinite(phi_var) || phi_var <= 0) {
      Rcpp::stop("Invalid phi_var at j = %d", j);
    }
    
    // ---- Posterior mean ----
    double phi_mean = phi_var * (
      (1.0 / sigma2[j]) * sum_FF_resid +
        mu_phi / sigma_phi2
    );
    
    double sd = std::sqrt(phi_var);
    
    // ---- Sample ----
    phi[free_idx] = rtnorm_lower_single_cpp(phi_mean, sd, 0);
    
    free_idx++;
  }
  
  // ---- Final consistency check ----
  if (free_idx != phi_factor.size()) {
    Rcpp::stop("phi_factor length does not match number of free phi");
  }
  
  return phi;
}

// [[Rcpp::export]]
Eigen::MatrixXd sample_Psi_cpp(
    const Eigen::VectorXd& beta,
    const Eigen::MatrixXd& n_vk,
    const IntegerVector& anchor_words, // 1-based
    double anchor_prob,
    double nonanchor_prob
) {
  
  int V = n_vk.rows();
  int K = n_vk.cols();
  
  double scale = 1.0 - anchor_prob - (K - 1) * nonanchor_prob;
  
  // ---- Mark anchor words ----
  std::vector<bool> is_anchor(V, false);
  
  for (int i = 0; i < anchor_words.size(); ++i) {
    int idx = anchor_words[i] - 1;
    if (idx < 0 || idx >= V) {
      Rcpp::stop("anchor_words index out of bounds");
    }
    is_anchor[idx] = true;
  }
  
  int VmK = V - anchor_words.size();
  
  Eigen::MatrixXd Psi(VmK, K);
  
  // ---- Build index map for free words ----
  std::vector<int> free_idx;
  free_idx.reserve(VmK);
  
  for (int v = 0; v < V; ++v) {
    if (!is_anchor[v]) {
      free_idx.push_back(v);
    }
  }
  
  // ---- Sample each column ----
  for (int k = 0; k < K; ++k) {
    
    Eigen::VectorXd gamma_draw(VmK);
    
    // ---- Draw gamma ----
    for (int i = 0; i < VmK; ++i) {
      int v = free_idx[i];
      double shape = beta[i] + n_vk(v, k);
      
      if (shape <= 0) {
        Rcpp::stop("Invalid Dirichlet parameter");
      }
      
      gamma_draw[i] = R::rgamma(shape, 1.0);
    }
    
    // ---- Normalize ----
    double sum_gamma = gamma_draw.sum();
    gamma_draw /= sum_gamma;
    
    // ---- Apply scaling ----
    Psi.col(k) = scale * gamma_draw;
  }
  
  return Psi;
}

// [[Rcpp::export]]
List sample_Z_cpp(
    List W,
    const NumericMatrix& Theta,
    const NumericMatrix& log_Psi_full   // <-- pass log(Psi_full)
) {
  
  int N = W.size();
  int K = Theta.ncol();
  
  List Z(N);
  
  // ---- reusable buffers ----
  std::vector<double> log_theta_i(K);
  std::vector<double> log_prob(K);
  std::vector<double> probs(K);
  
  for (int i = 0; i < N; i++) {
    
    IntegerVector Wi = W[i];
    int Mi = Wi.size();
    
    IntegerVector Zi(Mi);
    
    // ---- precompute log(theta_i) ----
    for (int k = 0; k < K; k++) {
      log_theta_i[k] = std::log(Theta(i, k));
    }
    
    for (int m = 0; m < Mi; m++) {
      
      int word = Wi[m] - 1;  // 0-based
      
      // ---- compute log probs + max in one pass ----
      double max_log = -INFINITY;
      
      for (int k = 0; k < K; k++) {
        double val = log_theta_i[k] + log_Psi_full(word, k);
        log_prob[k] = val;
        if (val > max_log) max_log = val;
      }
      
      // ---- exponentiate + sum ----
      double sum_probs = 0.0;
      
      for (int k = 0; k < K; k++) {
        double p = std::exp(log_prob[k] - max_log);
        probs[k] = p;
        sum_probs += p;
      }
      
      // ---- sample categorical ----
      double u = R::runif(0.0, sum_probs);
      double cum = 0.0;
      
      int chosen = K - 1;  // fallback (numerical safety)
      
      for (int k = 0; k < K; k++) {
        cum += probs[k];
        if (u <= cum) {
          chosen = k;
          break;
        }
      }
      
      Zi[m] = chosen + 1;  // back to 1-based
    }
    
    Z[i] = Zi;
  }
  
  return Z;
}

// [[Rcpp::export]]
Eigen::MatrixXd comp_Theta_cpp(
    const Eigen::VectorXd& tau_full, // Length K     
    const Eigen::MatrixXd& FF,     // N x L
    const Eigen::MatrixXd& Gamma  // L x K
) {
  
  int N = FF.rows();
  int K = Gamma.cols();
  
  Eigen::MatrixXd Theta(N, K);
  
  // ---- Compute linear predictor ----
  Eigen::MatrixXd lin = FF * Gamma;  // N x K
  
  // ---- Add Tau (broadcast across rows) ----
  lin.rowwise() += tau_full.transpose();
  
  // ---- Row-wise softmax with stability ----
  for (int i = 0; i < N; ++i) {
    
    // find max
    double max_val = lin(i, 0);
    for (int k = 1; k < K; ++k) {
      if (lin(i, k) > max_val) max_val = lin(i, k);
    }
    
    // compute exp + sum
    double sum_exp = 0.0;
    
    for (int k = 0; k < K; ++k) {
      double val = std::exp(lin(i, k) - max_val);
      Theta(i, k) = val;
      sum_exp += val;
    }
    
    // normalize
    for (int k = 0; k < K; ++k) {
      Theta(i, k) /= sum_exp;
    }
  }
  
  return Theta;
}

// [[Rcpp::export]]
NumericMatrix comp_n_vk_cpp(
    List W,
    List Z,
    int V,
    int K
) {
  
  int N = W.size();
  NumericMatrix n_vk(V, K);
  
  for (int i = 0; i < N; ++i) {
    
    IntegerVector w_i = W[i];
    IntegerVector z_i = Z[i];
    int M_i = w_i.size();
    
    for (int m = 0; m < M_i; ++m) {
      
      int v = w_i[m] - 1;  // convert to 0-based
      int k = z_i[m] - 1;
      
      n_vk(v, k) += 1.0;
    }
  }
  
  return n_vk;
}

// [[Rcpp::export]]
List make_Z_cpp(NumericMatrix Theta, IntegerVector M) {
  
  int N = Theta.nrow();
  int K = Theta.ncol();
  
  List Z(N);
  
  for (int i = 0; i < N; ++i) {
    
    // build CDF
    NumericVector cs(K);
    cs[0] = Theta(i, 0);
    for (int k = 1; k < K; ++k) {
      cs[k] = cs[k-1] + Theta(i, k);
    }
    
    int Mi = M[i];
    IntegerVector draws(Mi);
    
    for (int m = 0; m < Mi; ++m) {
      double u = R::runif(0.0, 1.0);
      
      int k = 0;
      while (u > cs[k]) ++k;
      
      draws[m] = k + 1;
    }
    
    Z[i] = draws;
  }
  
  return Z;
}

// [[Rcpp::export]]
List make_W_cpp(const NumericMatrix& Psi_full, const List& Z) {
  int N = Z.size();
  int V = Psi_full.nrow();
  
  List W(N);
  
  for (int i = 0; i < N; i++) {
    IntegerVector Zi = Z[i];
    int Mi = Zi.size();
    
    NumericVector Wi(Mi);
    
    for (int m = 0; m < Mi; m++) {
      int z = Zi[m] - 1;  // convert to 0-based index
      
      // Draw from categorical distribution
      double u = R::runif(0.0, 1.0);
      double cum_prob = 0.0;
      
      for (int v = 0; v < V; v++) {
        cum_prob += Psi_full(v, z);
        if (u <= cum_prob) {
          Wi[m] = v + 1;  // back to 1-based index
          break;
        }
      }
    }
    
    W[i] = Wi;
  }
  
  return W;
}

// [[Rcpp::export]]
Eigen::VectorXd comp_log_p_FF_cpp(
    const Eigen::MatrixXd& Y,        // N x J
    const Eigen::MatrixXd& n_ik,     // N x K
    const Eigen::MatrixXd& FF,       // N x L
    const Eigen::VectorXd& mu,       // length J
    const Eigen::MatrixXd& R_inv,    // L x L
    const Eigen::VectorXd& sigma2,   // length J  <-- updated
    const Eigen::MatrixXd& Phi,      // L x J
    const Eigen::MatrixXd& Gamma,    // L x K
    const Eigen::VectorXd& tau,      // length K
    const Eigen::MatrixXd& FF_Gamma  // N x K
) {
  
  int N = FF.rows();
  int J = Y.cols();
  int K = Gamma.cols();
  
  Eigen::VectorXd log_p(N);
  
  // ---- Precompute ----
  Eigen::MatrixXd FF_Rinv = FF * R_inv;          // N x L
  Eigen::MatrixXd FF_Phi  = FF * Phi;            // N x J
  Eigen::MatrixXd nG = n_ik * Gamma.transpose(); // N x L
  Eigen::VectorXd M = n_ik.rowwise().sum();      // N
  
  // ---- Main loop ----
  for (int i = 0; i < N; ++i) {
    
    // ---- Prior ----
    double log_prior = -0.5 * FF.row(i).dot(FF_Rinv.row(i));
    
    // ---- CFA likelihood ----
    double log_lik_cfa = 0.0;
    for (int j = 0; j < J; ++j) {
      double resid = Y(i, j) - mu[j] - FF_Phi(i, j);
      log_lik_cfa += (resid * resid) / sigma2[j];
    }
    log_lik_cfa *= -0.5;
    
    // ---- LDA numerator ----
    double log_lik_num = FF.row(i).dot(nG.row(i));
    
    // ---- LDA denominator (log-sum-exp) ----
    double max_val = FF_Gamma(i, 0) + tau[0];
    
    for (int k = 1; k < K; ++k) {
      double val = FF_Gamma(i, k) + tau[k];
      if (val > max_val) max_val = val;
    }
    
    double sum_exp = 0.0;
    for (int k = 0; k < K; ++k) {
      double val = FF_Gamma(i, k) + tau[k];
      sum_exp += std::exp(val - max_val);
    }
    
    double log_lik_den = M[i] * (max_val + std::log(sum_exp));
    
    // ---- Combine ----
    log_p[i] = log_prior + log_lik_cfa + log_lik_num - log_lik_den;
  }
  
  return log_p;
}

// [[Rcpp::export]]
Eigen::MatrixXd comp_grad_log_p_FF_cpp(
    const Eigen::MatrixXd& Y,        // N x J
    const Eigen::MatrixXd& n_ik,     // N x K
    const Eigen::MatrixXd& FF,       // N x L
    const Eigen::VectorXd& mu,       // length J
    const Eigen::MatrixXd& R_inv,    // L x L
    const Eigen::VectorXd& sigma2,   // length J
    const Eigen::MatrixXd& Phi,      // L x J
    const Eigen::MatrixXd& Gamma,    // L x K
    const Eigen::VectorXd& tau_full, // length K
    const Eigen::MatrixXd& FF_Gamma  // N x K
) {
  
  int N = FF.rows();
  int J = Y.cols();
  int K = Gamma.cols();
  int L = FF.cols();
  
  Eigen::MatrixXd grad(N, L);
  
  // ---- Precompute ----
  Eigen::MatrixXd grad_prior = -FF * R_inv;     // N x L
  Eigen::VectorXd M = n_ik.rowwise().sum();    // N
  
  // ---- Main loop ----
  for (int i = 0; i < N; ++i) {
    
    // ---- CFA gradient ----
    Eigen::VectorXd grad_cfa_i = Eigen::VectorXd::Zero(L);
    
    for (int j = 0; j < J; ++j) {
      double resid = Y(i, j) - mu[j];
      
      // subtract FF * Phi contribution
      for (int l = 0; l < L; ++l) {
        resid -= FF(i, l) * Phi(l, j);
      }
      
      double scaled = resid / sigma2[j];
      
      for (int l = 0; l < L; ++l) {
        grad_cfa_i[l] += scaled * Phi(l, j);
      }
    }
    
    // ---- Softmax (Theta row) ----
    double max_val = FF_Gamma(i, 0) + tau_full[0];
    for (int k = 1; k < K; ++k) {
      double val = FF_Gamma(i, k) + tau_full[k];
      if (val > max_val) max_val = val;
    }
    
    Eigen::VectorXd theta_i(K);
    double sum_exp = 0.0;
    
    for (int k = 0; k < K; ++k) {
      double val = std::exp(FF_Gamma(i, k) + tau_full[k] - max_val);
      theta_i[k] = val;
      sum_exp += val;
    }
    
    theta_i /= sum_exp;
    
    // ---- LDA gradient ----
    Eigen::VectorXd lda_term = n_ik.row(i).transpose() - M[i] * theta_i;
    Eigen::VectorXd grad_lda_i = Gamma * lda_term;
    
    // ---- Combine ----
    grad.row(i) = grad_prior.row(i) + grad_cfa_i.transpose() + grad_lda_i.transpose();
  }
  
  return grad;
}

// [[Rcpp::export]]
double comp_log_p_gamma_cpp(
    const Eigen::MatrixXd& FF_Gamma, // N x K (precomputed)
    const Eigen::VectorXd& gamma,    // vectorized Gamma
    const Eigen::MatrixXd& n_ik,     // N x K
    const Eigen::VectorXd& tau_full, // length K
    double mu_gamma,                 // scalar
    double sigma_gamma2
) {
  
  int N = FF_Gamma.rows();
  int K = FF_Gamma.cols();
  
  Eigen::VectorXd M = n_ik.rowwise().sum();
  
  double loglik_num = 0.0;
  double loglik_den = 0.0;
  
  // ---- Main loop ----
  for (int i = 0; i < N; ++i) {
    
    double max_val = FF_Gamma(i, 0) + tau_full[0];
    
    for (int k = 1; k < K; ++k) {
      double val = FF_Gamma(i, k) + tau_full[k];
      if (val > max_val) max_val = val;
    }
    
    double sum_exp = 0.0;
    double row_num = 0.0;
    
    for (int k = 0; k < K; ++k) {
      double val = FF_Gamma(i, k) + tau_full[k];
      
      row_num += n_ik(i, k) * val;
      sum_exp += std::exp(val - max_val);
    }
    
    loglik_num += row_num;
    loglik_den += M[i] * (max_val + std::log(sum_exp));
  }
  
  // ---- Prior ----
  double sq_sum = 0.0;
  for (int i = 0; i < gamma.size(); ++i) {
    double diff = gamma[i] - mu_gamma;
    sq_sum += diff * diff;
  }
  
  double log_prior = -0.5 * sigma_gamma2 * sq_sum;
  
  return loglik_num - loglik_den + log_prior;
}

// [[Rcpp::export]]
Eigen::VectorXd comp_grad_log_p_gamma_cpp(
    const Eigen::MatrixXd& FF,        // N x L
    const Eigen::MatrixXd& FF_Gamma,  // N x K
    const Eigen::VectorXd& gamma,     // vectorized Gamma (subset)
    const Eigen::MatrixXd& n_ik,      // N x K
    const Eigen::VectorXd& tau_full,       // length K
    double mu_gamma,                  // scalar
    double sigma_gamma2,
    const IntegerVector& gamma_fixed, // 1-based (rows of Gamma)
    const IntegerVector& gamma_factor // 1-based (cols of Gamma)
) {
  
  int N = FF.rows();
  int L = FF.cols();
  int K = FF_Gamma.cols();
  int Pdim = gamma.size();
  
  Eigen::MatrixXd grad_Gamma = Eigen::MatrixXd::Zero(L, K);
  
  Eigen::VectorXd M = n_ik.rowwise().sum();
  
  // ---- Main loop over observations ----
  for (int i = 0; i < N; ++i) {
    
    // ---- softmax ----
    double max_val = FF_Gamma(i, 0) + tau_full[0];
    for (int k = 1; k < K; ++k) {
      double val = FF_Gamma(i, k) + tau_full[k];
      if (val > max_val) max_val = val;
    }
    
    Eigen::VectorXd P(K);
    double sum_exp = 0.0;
    
    for (int k = 0; k < K; ++k) {
      double val = std::exp(FF_Gamma(i, k) + tau_full[k] - max_val);
      P[k] = val;
      sum_exp += val;
    }
    
    P /= sum_exp;
    
    // ---- compute (n_ik - M * P) ----
    Eigen::VectorXd diff = n_ik.row(i).transpose() - M[i] * P;
    
    // ---- accumulate gradient: FF_i^T * diff ----
    // outer product: (L x 1) * (1 x K)
    grad_Gamma += FF.row(i).transpose() * diff.transpose();
  }
  
  // ---- Extract gradient for gamma subset ----
  Eigen::VectorXd grad_gamma(Pdim);
  
  for (int p = 0; p < Pdim; ++p) {
    int l = gamma_fixed[p] - 1;   // row index
    int k = gamma_factor[p] - 1;  // column index
    
    if (l < 0 || l >= L || k < 0 || k >= K) {
      Rcpp::stop("gamma indices out of bounds");
    }
    
    grad_gamma[p] = grad_Gamma(l, k);
  }
  
  // ---- Prior ----
  double sum_diff = 0.0;
  for (int i = 0; i < Pdim; ++i) {
    sum_diff += (gamma[i] - mu_gamma);
  }
  
  grad_gamma.array() += sum_diff / sigma_gamma2;
  
  return grad_gamma;
}

// [[Rcpp::export]]
double comp_log_p_tau_cpp(
    const Eigen::MatrixXd& FF_Gamma, // N x K
    const Eigen::MatrixXd& n_ik,     // N x K
    const Eigen::VectorXd& tau,      // length K-1
    const Eigen::VectorXd& tau_full, // length K
    double mu_tau,                   // scalar
    double sigma_tau2
) {
  
  int N = FF_Gamma.rows();
  int K = FF_Gamma.cols();
  
  Eigen::VectorXd M = n_ik.rowwise().sum();
  
  double loglik_num = 0.0;
  double loglik_den = 0.0;
  
  // ---- Main loop ----
  for (int i = 0; i < N; ++i) {
    
    // ---- log-sum-exp max ----
    double max_val = FF_Gamma(i, 0) + tau_full[0];
    
    for (int k = 1; k < K; ++k) {
      double val = FF_Gamma(i, k) + tau_full[k];
      if (val > max_val) max_val = val;
    }
    
    double sum_exp = 0.0;
    double row_num = 0.0;
    
    for (int k = 0; k < K; ++k) {
      double val = FF_Gamma(i, k) + tau_full[k];
      
      row_num += n_ik(i, k) * val;
      sum_exp += std::exp(val - max_val);
    }
    
    loglik_num += row_num;
    loglik_den += M[i] * (max_val + std::log(sum_exp));
  }
  
  // ---- Prior ----
  double sq_sum = 0.0;
  for (int k = 0; k < tau.size(); ++k) {
    double diff = tau[k] - mu_tau;
    sq_sum += diff * diff;
  }
  
  double log_prior = -0.5 * sigma_tau2 * sq_sum;
  
  return loglik_num - loglik_den + log_prior;
}

// [[Rcpp::export]]
Eigen::VectorXd comp_grad_log_p_tau_cpp(
    const Eigen::MatrixXd& FF_Gamma, // N x K
    const Eigen::MatrixXd& n_ik,     // N x K
    const Eigen::VectorXd& tau,      // length K-1
    const Eigen::VectorXd& tau_full, // length K
    double mu_tau,                   // scalar
    double sigma_tau2,
    const IntegerVector& topic_baseline // 1-based index to drop
) {
  
  int N = FF_Gamma.rows();
  int K = FF_Gamma.cols();
  
  Eigen::VectorXd M = n_ik.rowwise().sum();   // N
  Eigen::VectorXd n = n_ik.colwise().sum();   // K
  
  Eigen::VectorXd grad_full = Eigen::VectorXd::Zero(K);
  
  // ---- Main loop ----
  for (int i = 0; i < N; ++i) {
    
    // ---- log-sum-exp ----
    double max_val = FF_Gamma(i, 0) + tau_full[0];
    
    for (int k = 1; k < K; ++k) {
      double val = FF_Gamma(i, k) + tau_full[k];
      if (val > max_val) max_val = val;
    }
    
    double sum_exp = 0.0;
    std::vector<double> exp_vals(K);
    
    for (int k = 0; k < K; ++k) {
      double val = std::exp(FF_Gamma(i, k) + tau_full[k] - max_val);
      exp_vals[k] = val;
      sum_exp += val;
    }
    
    // ---- accumulate gradient ----
    for (int k = 0; k < K; ++k) {
      double P_ik = exp_vals[k] / sum_exp;
      grad_full[k] -= M[i] * P_ik;
    }
  }
  
  // ---- add n ----
  grad_full += n;
  
  // ---- Prior ----
  double sum_diff = 0.0;
  for (int k = 0; k < tau.size(); ++k) {
    sum_diff += (tau[k] - mu_tau);
  }
  
  grad_full.array() += sum_diff / sigma_tau2;
  
  // ---- Remove baseline ----
  int b = topic_baseline[0] - 1;
  
  Eigen::VectorXd grad(K - 1);
  int idx = 0;
  
  for (int k = 0; k < K; ++k) {
    if (k == b) continue;
    grad[idx++] = grad_full[k];
  }
  
  return grad;
}