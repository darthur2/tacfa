#include <RcppEigen.h>

// Rcpp::depends(Rcpp)
// Rcpp::depends(RcppEigen)

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
double comp_ldelta_k_cpp(double delta_k,
                         const Eigen::MatrixXd& FF,
                         const Eigen::MatrixXd& Gamma,
                         const Eigen::MatrixXd& Z_sum,
                         const Eigen::VectorXd& delta,
                         int k){
  
  Eigen::VectorXd delta_k_vec = delta;
  delta_k_vec(k-1) = delta_k;
  
  Eigen::MatrixXd FF_Gamma_k = (FF*Gamma).rowwise() + delta_k_vec.transpose();
  Eigen::VectorXd FFG_log_exp_sums = FF_Gamma_k.array().exp().rowwise().sum().log();
  
  double ldelta_k = (Z_sum.array()*((FF_Gamma_k.colwise() - FFG_log_exp_sums).array())).sum() - 
    0.5*pow(delta_k, 2);
  
  return ldelta_k;
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