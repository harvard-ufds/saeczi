// [[Rcpp::depends(RcppEigen)]]

#include <string>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "predict.h"
using namespace Eigen;
using namespace Rcpp;


//[[Rcpp::export]]
SEXP predict(const Rcpp::DataFrame &newdata,
             const String domain,
             const Eigen::MatrixXd &lm,
             const Eigen::MatrixXd &glm,
             const Rcpp::List &u,
             const Rcpp::CharacterVector lm_X,
             const Rcpp::CharacterVector glm_X,
             const int B) {
  
  
  
  // mat_u_lin <- matrix(rep(0, N*B), nrow = N)
  // mat_u_log <- matrix(rep(0, N*B), nrow = N)
  int N = newdata.nrows();
  Eigen::MatrixXd mat_u_lm(N, B);
  Eigen::MatrixXd mat_u_glm(N, B);
  
  for (int j = 0; j < B; ++j) {
    Eigen::VectorXd colVec_lm = as<Eigen::VectorXd>(u[j]);
    Eigen::VectorXd colVec_glm = as<Eigen::VectorXd>(u[j]);
    mat_u_lm.col(j) = colVec_lm;
    mat_u_glm.col(j) = colVec_glm;
    
  }
  // for (i in seq_len(length(u))) {
  //  mat_u_lin[ ,i] <- u[[i]]$u_lm[dom_ref]
  //  mat_u_log[ ,i] <- u[[i]]$u_glm[dom_ref]
  //}
  
  
  return(0);
}


