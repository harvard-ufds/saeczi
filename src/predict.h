#ifndef predict_H
#define predict_H

#include <string>
#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Eigen;
using namespace Rcpp;

SEXP predict(const Rcpp::DataFrame &newdata,
             const String &domain,
             const Eigen::MatrixXd &lm,
             const Eigen::MatrixXd &glm,
             const Rcpp::List &u,
             const Rcpp::StringVector &lm_X,
             const Rcpp::StringVector &glm_X);

#endif