#include <RcppEigen.h>

// Rcpp::depends(Rcpp)
// Rcpp::depends(RcppEigen)

// [[Rcpp::export]]
Eigen::MatrixXd sum_Z_cpp(const std::vector<Eigen::VectorXi>& Z,
                          int K) {
  int D = Z.size();
  
  Eigen::MatrixXd Z_sum = Eigen::MatrixXd::Zero(D, K);
  
  for (int d = 0; d < D; ++d){
    Eigen::VectorXi z_d = Z[d];
    
    for (int i = 0; i < z_d.size(); ++i){
      int k = z_d(i)-1;
      
      Z_sum(d,k) += 1;
    }
  }
  
  return Z_sum;
}

// [[Rcpp::export]]
Eigen::MatrixXd rdirichlet_cpp(int n,
                               const Eigen::VectorXd& shape){
  int m = shape.size();
  
  Eigen::MatrixXd draws(m, n);
  
  for (int i = 0; i < n; ++i){
    Eigen::VectorXd gamma_draws_i(m);
    
    for (int j = 0; j < m; ++j){
      gamma_draws_i(j) = R::rgamma(shape[j], 1.0);
    }
    
    draws.col(i) = gamma_draws_i / gamma_draws_i.sum();
  }
  
  return draws;
}

// [[Rcpp::export]]
Eigen::MatrixXd comp_WZ_counts_cpp(const std::vector<Eigen::VectorXi>& W,
                                   const std::vector<Eigen::VectorXi>& Z,
                                   int V,
                                   int K){
  
  int D = W.size();
  Eigen::MatrixXd WZ_counts = Eigen::MatrixXd::Zero(V, K);
  
  for (int d = 0; d < D; ++d){
    Eigen::VectorXi w_d = W[d];
    Eigen::VectorXi z_d = Z[d];
    
    for (int i = 0; i < w_d.size(); ++i){
      int word = w_d[i]-1;
      int topic = z_d[i]-1;
      
      WZ_counts(word, topic) += 1;
    }
  }
  
  return WZ_counts;
}

// [[Rcpp::export]]
Eigen::MatrixXd rmvnorm_cpp(int n,
                            const Eigen::VectorXd& Mu,
                            const Eigen::MatrixXd& Sigma){
  
  // Cholesky decomposition
  Eigen::LLT<Eigen::MatrixXd> lltOfSigma(Sigma); // compute the Cholesky decomposition of Sigma
  Eigen::MatrixXd Q = lltOfSigma.matrixL().transpose(); // get the lower triangular matrix L such that Sigma = L * L^T
  int L = Mu.size();
  
  Eigen::MatrixXd X = Eigen::MatrixXd::NullaryExpr(n, L, [&](){ return R::rnorm(0.0, 1.0); });
  Eigen::MatrixXd Y = (X * Q).rowwise() + Mu.transpose();
  
  return Y;
}

// [[Rcpp::export]]
Eigen::VectorXd rtnorm_cpp(int n,
                           double mu,
                           double sigma,
                           double lower,
                           double upper){
  
  double p_lower = R::pnorm(lower, mu, sigma, 1, 0);
  double p_upper = R::pnorm(upper, mu, sigma, 1, 0);
  
  Eigen::VectorXd u = Eigen::VectorXd::NullaryExpr(n, [=](){ 
    return R::runif(p_lower, p_upper); 
    });
  
  Eigen::VectorXd draws = u.unaryExpr([=](double ui) {
    return R::qnorm5(ui, mu, sigma, 1, 0);
  });
  
  return draws;
}

// [[Rcpp::export]]
int sample_cpp(const Eigen::VectorXd& p) {
  double u = R::runif(0.0, 1.0);
  double cumulative = 0.0;
  
  for (int i = 0; i < p.size(); ++i) {
    cumulative += p[i];
    if (u <= cumulative) {
      return i;
    }
  }
  
  return p.size();  // fallback for numerical error
}

// [[Rcpp::export]]
double comp_lgamma_lk_cpp(double gamma_lk,
                          const Eigen::MatrixXd& FF,
                          const Eigen::MatrixXd& Gamma,
                          const Eigen::MatrixXd& Z_sum,
                          const Eigen::VectorXd& delta,
                          int l,
                          int k){
  
  
  Eigen::MatrixXd Gamma_lk = Gamma;
  Gamma_lk(l,k) = gamma_lk;
  
  Eigen::MatrixXd FF_Gamma_lk = (FF*Gamma_lk).rowwise() + delta.transpose();
  Eigen::VectorXd FFG_log_exp_sums = FF_Gamma_lk.array().exp().rowwise().sum().log();
  
  double lgamma_lk = (Z_sum.array()*((FF_Gamma_lk.colwise() - FFG_log_exp_sums).array())).sum() - 
    0.5*pow(gamma_lk, 2);
  
  return lgamma_lk;
}

// [[Rcpp::export]]
double sample_gamma_lk_cpp(const Eigen::MatrixXd& FF,
                           const Eigen::MatrixXd& Gamma,
                           const Eigen::MatrixXd& Z_sum,
                           const Eigen::VectorXd& delta,
                           double gamma_sd,
                           int l,
                           int k) {
  
  double gamma_lk = Gamma(l, k);
  Eigen::MatrixXd Gamma_star = Gamma;
  double gamma_lk_star = R::rnorm(gamma_lk, gamma_sd);  // Propose new value
  Gamma_star(l, k) = gamma_lk_star;
  
  double log_gamma_lk = comp_lgamma_lk_cpp(gamma_lk, FF, Gamma, Z_sum, delta, l, k);
  double log_gamma_lk_star = comp_lgamma_lk_cpp(gamma_lk_star, FF, Gamma_star, Z_sum, delta, l, k);
  
  double log_u = std::log(R::runif(0.0, 1.0));
  double log_alpha = log_gamma_lk_star - log_gamma_lk;
  
  if (log_alpha > log_u) {
    return gamma_lk_star;
  } else {
    return gamma_lk;
  }
}

// [[Rcpp::export]]
Eigen::MatrixXd sample_Gamma_cpp(const Eigen::MatrixXd& FF,
                                 Eigen::MatrixXd Gamma,
                                 const Eigen::MatrixXd& Z_sum,
                                 const Eigen::VectorXd& delta,
                                 double gamma_sd) {
  
  int L = Gamma.rows();
  int K = Gamma.cols();
  
  for (int l = 0; l < L; ++l) {
    for (int k = 0; k < K - 1; ++k) {  // K-1 due to constraint (e.g., identifiability)
      Gamma(l, k) = sample_gamma_lk_cpp(FF, Gamma, Z_sum, delta, gamma_sd, l, k);
    }
  }
  
  return Gamma;
}

// [[Rcpp::export]]
double comp_lf_d_cpp(const Eigen::VectorXd& f_d,
                     const Eigen::VectorXd& Y_d,
                     const Eigen::VectorXd& Z_sum_d,
                     const Eigen::MatrixXd& Phi,
                     const Eigen::MatrixXd& Gamma,
                     const Eigen::MatrixXd& R_inv,
                     const Eigen::VectorXd& mu,
                     const Eigen::VectorXd& sigma2,
                     const Eigen::VectorXd& delta){
  
  
  Eigen::MatrixXd Phi_sigma = Phi.array().rowwise() / sigma2.array().sqrt().transpose();
  
  Eigen::MatrixXd f_d_Sigma_inv = Phi_sigma*Phi_sigma.transpose() + R_inv;
  
  Eigen::VectorXd weights = (Y_d - mu).array() / sigma2.array();;
  Eigen::VectorXd f_d_Mu = Phi * weights;
  
  Eigen::VectorXd f_d_Gamma = (f_d.transpose() * Gamma + delta.transpose()).transpose();
  Eigen::VectorXd exp_f_d_Gamma = f_d_Gamma.array().exp();
  Eigen::VectorXd theta_d = exp_f_d_Gamma / exp_f_d_Gamma.sum();
  
  double Z_sum_log_theta_d = (Z_sum_d.array()*theta_d.array().log()).sum();
  double quad_term = (f_d.transpose() * f_d_Sigma_inv * f_d)(0, 0);
  double linear_term = (f_d.transpose() * f_d_Mu)(0, 0);
  double lf_d = -0.5 * (quad_term - 2.0 * linear_term) + Z_sum_log_theta_d;
  
  return lf_d;
}

// [[Rcpp::export]]
Eigen::VectorXd sample_f_d_cpp(const Eigen::VectorXd& f_d,
                               const Eigen::VectorXd& Y_d,
                               const Eigen::VectorXd& Z_sum_d,
                               const Eigen::MatrixXd& Phi,
                               const Eigen::MatrixXd& Gamma,
                               const Eigen::MatrixXd& R_inv,
                               const Eigen::VectorXd& mu,
                               const Eigen::VectorXd& sigma2,
                               const Eigen::VectorXd& delta,
                               const Eigen::MatrixXd& f_Sigma){
  
  Eigen::VectorXd f_d_star = rmvnorm_cpp(1, f_d, f_Sigma).row(0);
  
  double log_f_d = comp_lf_d_cpp(f_d, Y_d, Z_sum_d, Phi, Gamma, 
                                 R_inv, mu, sigma2, delta);
  double log_f_d_star = comp_lf_d_cpp(f_d_star, Y_d, Z_sum_d, Phi, Gamma, 
                                      R_inv, mu, sigma2, delta);
  
  double log_u = std::log(R::runif(0.0, 1.0));
  double log_alpha = log_f_d_star - log_f_d;
  
  if (log_alpha > log_u){
    return f_d_star;
  } else {
    return f_d;
  }
}

// [[Rcpp::export]]
Eigen::MatrixXd sample_FF_cpp(Eigen::MatrixXd FF,
                              const Eigen::MatrixXd& Y,
                              const Eigen::MatrixXd& Z_sum,
                              const Eigen::MatrixXd& Phi,
                              const Eigen::MatrixXd& Gamma,
                              const Eigen::MatrixXd& R_inv,
                              const Eigen::VectorXd& mu,
                              const Eigen::VectorXd& sigma2,
                              const Eigen::VectorXd& delta,
                              const Eigen::MatrixXd f_Sigma){
  
  int D = FF.rows();
  
  for (int d = 0; d < D; ++d){
    FF.row(d) = sample_f_d_cpp(FF.row(d), Y.row(d), Z_sum.row(d), Phi, Gamma,
           R_inv, mu, sigma2, delta, f_Sigma);
  }
  
  return FF;
}

// [[Rcpp::export]]
double comp_ldelta_k_cpp(double delta_k,
                         const Eigen::MatrixXd& FF,
                         const Eigen::MatrixXd& Gamma,
                         const Eigen::MatrixXd& Z_sum,
                         const Eigen::VectorXd& delta,
                         int k){
  
  Eigen::VectorXd delta_k_vec = delta;
  delta_k_vec(k) = delta_k;
  
  Eigen::MatrixXd FF_Gamma_k = (FF*Gamma).rowwise() + delta_k_vec.transpose();
  Eigen::VectorXd FFG_log_exp_sums = FF_Gamma_k.array().exp().rowwise().sum().log();
  
  double ldelta_k = (Z_sum.array()*((FF_Gamma_k.colwise() - FFG_log_exp_sums).array())).sum() - 
    0.5*pow(delta_k, 2);
  
  return ldelta_k;
}

// [[Rcpp::export]]
double sample_delta_k_cpp(const Eigen::MatrixXd& FF,
                          const Eigen::MatrixXd& Gamma,
                          const Eigen::MatrixXd& Z_sum,
                          const Eigen::VectorXd& delta,
                          double delta_sd,
                          int k){
  double delta_k = delta(k);
  double delta_k_star = R::rnorm(delta_k, delta_sd);
  
  double log_delta_k = comp_ldelta_k_cpp(delta_k, FF, Gamma, Z_sum, delta, k);
  double log_delta_k_star = comp_ldelta_k_cpp(delta_k_star, FF, Gamma, Z_sum, delta, k);
  
  double log_u = std::log(R::runif(0.0, 1.0));
  double log_alpha = log_delta_k_star - log_delta_k;
  
  if (log_alpha > log_u){
    return delta_k_star;
  } else {
    return delta_k;
  }
}

// [[Rcpp::export]]
Eigen::VectorXd sample_delta_cpp(const Eigen::MatrixXd& FF,
                                 const Eigen::MatrixXd& Gamma,
                                 const Eigen::MatrixXd& Z_sum,
                                 Eigen::VectorXd delta,
                                 double delta_sd){
  int K = delta.size();
  
  for (int k = 0; k < K-1; ++k){
    delta(k) = sample_delta_k_cpp(FF, Gamma, Z_sum, delta, delta_sd, k);
  }
  
  return delta;
}

// [[Rcpp::export]]
std::vector<Eigen::MatrixXd> comp_Z_post_cpp(const std::vector<Eigen::VectorXi>& W,
                                             const Eigen::MatrixXd& FF,
                                             const Eigen::MatrixXd& Gamma,
                                             const Eigen::MatrixXd& Psi,
                                             const Eigen::VectorXd& delta){
  int D = W.size();
  std::vector<Eigen::MatrixXd> Z_post(D);
  
  Eigen::MatrixXd FF_Gamma = FF * Gamma;
  FF_Gamma.rowwise() += delta.transpose();
  
  Eigen::MatrixXd exp_FF_Gamma = FF_Gamma.array().exp();
  Eigen::MatrixXd Theta = exp_FF_Gamma.array().colwise() / exp_FF_Gamma.rowwise().sum().array();
  
  for (int d = 0; d < D; ++d){
    int Md = W[d].size();
    int K = Psi.cols();
    Z_post[d] = Eigen::MatrixXd(Md, K);
    
    for (int m = 0; m < Md; ++m){
      int w_idx = W[d](m);
      Eigen::VectorXd psi_dm = Psi.row(w_idx-1);
      Eigen::VectorXd Z_post_unnorm_md = psi_dm.array()*Theta.row(d).transpose().array();
      Z_post[d].row(m) = Z_post_unnorm_md.transpose() / Z_post_unnorm_md.sum();
    }
  }
  
  return Z_post;
}

// [[Rcpp::export]]
std::vector<Eigen::VectorXi> sample_Z_cpp(std::vector<Eigen::VectorXi>& W,
                                          const Eigen::MatrixXd& FF,
                                          const Eigen::MatrixXd& Gamma,
                                          const Eigen::MatrixXd& Psi,
                                          const Eigen::VectorXd& delta){
  
  int D = W.size();
  int K = Gamma.cols();
  
  std::vector<Eigen::VectorXi> Z(D);
  std::vector<Eigen::MatrixXd> Z_post = comp_Z_post_cpp(W, FF, Gamma, Psi, delta);
  
  for (int d = 0; d < D; ++d){
    int Md = Z_post[d].rows();
    Z[d].resize(Md);
    
    for (int i = 0; i < Z_post[d].rows(); ++i){
      Z[d](i) = sample_cpp(Z_post[d].row(i)) + 1;
    }
  }
  
  return Z;
}

// [[Rcpp::export]]
double sample_phi_lj_cpp(const Eigen::MatrixXd& Y,
                         const Eigen::MatrixXd& FF,
                         const Eigen::VectorXd& mu,
                         const Eigen::VectorXd& sigma2,
                         const Eigen::MatrixXd& Phi,
                         int l,
                         int j){
  
  int D = FF.rows();
  int L = FF.cols();
  
  Eigen::MatrixXd FF_sub(D, L - 1);
  FF_sub << FF.leftCols(l), FF.rightCols(L - l - 1);  // skip col l
  
  Eigen::VectorXd Phi_sub(L - 1);
  Phi_sub << Phi.topRows(l).col(j), Phi.bottomRows(L - l - 1).col(j);  // skip row l
  
  Eigen::MatrixXd FF_Phi_sub = FF_sub * Phi_sub;
  
  double phi_lj_var = 1/(1/sigma2(j)*FF.col(l).array().square().sum() + 1.0);
  double phi_lj_mu = 1/sigma2(j)*((Y.col(j).array() - mu(j) - FF_Phi_sub.array())*FF.col(l).array()).sum();
  
  double phi_lj = rtnorm_cpp(1, phi_lj_mu*phi_lj_var, std::sqrt(phi_lj_var), 0, R_PosInf)(0);
  
  return phi_lj;
}

// [[Rcpp::export]]
Eigen::MatrixXd sample_Phi_cpp(const Eigen::MatrixXd& Y,
                               const Eigen::MatrixXd& FF,
                               const Eigen::VectorXd& mu,
                               const Eigen::VectorXd& sigma2,
                               Eigen::MatrixXd& Phi,
                               std::vector<Eigen::VectorXi>& Phi_indices){
  
  int L = FF.cols();
  
  for (int l = 0; l < L; ++l){
    for (int j = 0; j < Phi_indices[l].size(); ++j){
      int idx = Phi_indices[l](j) - 1;
      Phi(l,idx) = sample_phi_lj_cpp(Y, FF, mu, sigma2, Phi, l, idx);
    }
  }
  
  return Phi;
}

// [[Rcpp::export]]
double sample_sigma2_j_cpp(const Eigen::MatrixXd& Y,
                           const Eigen::VectorXd& mu,
                           const Eigen::MatrixXd& FF,
                           const Eigen::MatrixXd& Phi,
                           double a_sigma2,
                           double b_sigma2,
                           int j){
  
  int D = FF.rows();
  
  Eigen::VectorXd FF_Phi_j = (FF*Phi).col(j);
  
  double a_sigma2_star = D/2 + a_sigma2;
  double b_sigma2_star = 0.5*((Y.col(j).array() - mu(j) - FF_Phi_j.array()).square().sum()) + b_sigma2;
  
  double sigma2_j = 1.0/R::rgamma(a_sigma2_star, 1/b_sigma2_star);
  
  return sigma2_j;
}

// [[Rcpp::export]]
Eigen::VectorXd sample_sigma2_cpp(const Eigen::MatrixXd& Y,
                                  const Eigen::VectorXd& mu,
                                  const Eigen::MatrixXd& FF,
                                  const Eigen::MatrixXd& Phi,
                                  double a_sigma2,
                                  double b_sigma2){
  
  int J = Phi.cols();
  Eigen::VectorXd sigma2(J);
  
  for (int j = 0; j < J; ++j){
    sigma2(j) = sample_sigma2_j_cpp(Y, mu, FF, Phi, a_sigma2, b_sigma2, j);
  }
  
  return sigma2;
}

// [[Rcpp::export]]
double sample_mu_j_cpp(const Eigen::MatrixXd& Y,
                       const Eigen::VectorXd& sigma2,
                       const Eigen::MatrixXd& FF,
                       const Eigen::MatrixXd& Phi,
                       double sigma2_mu,
                       int j){
  
  int D = FF.rows();
  Eigen::VectorXd FF_Phi_j = (FF*Phi).col(j);
  
  double mu_j_var = 1/(D/sigma2(j) + 1/sigma2_mu);
  double mu_j_mu = 1/sigma2(j)*(Y.col(j).array() - FF_Phi_j.array()).sum();
  
  double mu_j = R::rnorm(mu_j_mu*mu_j_var, std::sqrt(mu_j_var));
  
  return mu_j;
}

// [[Rcpp::export]]
Eigen::VectorXd sample_mu_cpp(const Eigen::MatrixXd& Y,
                              const Eigen::VectorXd& sigma2,
                              const Eigen::MatrixXd& FF,
                              const Eigen::MatrixXd& Phi,
                              double sigma2_mu){
  
  int J = Phi.cols();
  Eigen::VectorXd mu(J);
  
  for (int j = 0; j < J; ++j){
    mu(j) = sample_mu_j_cpp(Y, sigma2, FF, Phi, sigma2_mu, j);
  }
  
  return mu;
}

// [[Rcpp::export]]
Eigen::MatrixXd sample_Psi_cpp(std::vector<Eigen::VectorXi>& W,
                               std::vector<Eigen::VectorXi>& Z,
                               const Eigen::MatrixXd& beta){
  
  int V = beta.rows();
  int K = beta.cols();
  
  Eigen::MatrixXd Psi(V, K);
  Eigen::MatrixXd WZ_counts = comp_WZ_counts_cpp(W, Z, V, K);
  Eigen::MatrixXd beta_star = WZ_counts + beta;
  
  for (int k = 0; k < K; ++k){
    Psi.col(k) = rdirichlet_cpp(1, beta_star.col(k));
  }
  
  return Psi;
}