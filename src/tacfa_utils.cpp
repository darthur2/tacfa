#include <RcppEigen.h>

// Rcpp::depends(Rcpp)
// Rcpp::depends(RcppEigen)

// [[Rcpp::export]]
Eigen::MatrixXd sum_Z_cpp(const std::vector<Eigen::VectorXi>& Z,
                          int K) {
  int D = Z.size();
  
  Eigen::MatrixXd Z_sum(D, K);
  
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
                               Eigen::VectorXd& shape){
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