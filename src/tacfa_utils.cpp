#include <RcppEigen.h>

// Rcpp::depends(Rcpp)
// Rcpp::depends(RcppEigen)

// [[Rcpp::export]]
Eigen::MatrixXd sum_Z_cpp(const std::vector<Eigen::SparseMatrix<double>>& Z) {
  int D = Z.size();
  int K = Z[0].cols();
  
  Eigen::MatrixXd Z_sum(D, K);
  
  for (int d = 0; d < D; ++d){
    Eigen::SparseMatrix<double> Zd = Z[d];
    Eigen::VectorXd Z_sum_d = Eigen::VectorXd::Zero(K);
    
    for (int k = 0; k < Zd.outerSize(); ++k) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(Zd, k); it; ++it) {
        Z_sum_d[it.col()] += it.value();
      }
    }
    
    Z_sum.row(d) = Z_sum_d;
  }
  
  return Z_sum;
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
  Gamma_lk(l-1,k-1) = gamma_lk;
  
  Eigen::MatrixXd FF_Gamma_lk = (FF*Gamma_lk).rowwise() + delta.transpose();
  Eigen::VectorXd FFG_log_exp_sums = FF_Gamma_lk.array().exp().rowwise().sum().log();
  
  double lgamma_lk = (Z_sum.array()*((FF_Gamma_lk.colwise() - FFG_log_exp_sums).array())).sum() - 
    0.5*pow(gamma_lk, 2);

  return lgamma_lk;
}

