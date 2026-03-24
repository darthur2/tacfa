#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int sample_catC(NumericVector p) {
  
  double u = R::runif(0.0, 1.0);
  double cum = 0.0;
  
  int K = p.size();
  
  for (int k = 0; k < K; k++) {
    cum += p[k];
    
    if (cum >= u) {
      return k + 1;  // R uses 1-based indexing
    }
  }
  
  // fallback (in case of numerical issues)
  return K;
}