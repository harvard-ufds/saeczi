// [[Rcpp::depends(RcppEigen)]]

#include <string>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "predict.h"
using namespace Eigen;
using namespace Rcpp;


//[[Rcpp::export]]
SEXP predict(const Rcpp::DataFrame &newdata,
             const String &domain,
             const Eigen::MatrixXd &lm,
             const Eigen::MatrixXd &glm,
             const Rcpp::List &u,
             const Rcpp::CharacterVector &lm_X,
             const Rcpp::CharacterVector &glm_X) {
  
  
  
  return(0);
}


