// [[Rcpp::depends(RcppEigen)]]
#include <string>
#include <unordered_map>
#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Eigen;
using namespace Rcpp;


//[[Rcpp::export]]
SEXP by_index(Rcpp::CharacterVector names,
              Rcpp::NumericVector vals,
              Rcpp::CharacterVector to_id) {
  
  std::unordered_map<Rcpp::String, double> namedVec;
  for(int i = 0; i < vals.size(); i ++) {
    namedVec[names[i]] = vals[i];
  }
  
  Rcpp::NumericVector res;
  
  for (const Rcpp::String& id : to_id) {
    res.push_back(namedVec[id]);
  }

  return Rcpp::wrap(res);
  
}

double sigmoid(double x) {
  return 1.0 / (1.0 + std::exp(-x));
}


//[[Rcpp::export]]
SEXP predict_zi(Eigen::MatrixXd &beta_lm,
                Eigen::MatrixXd &beta_glm,
                Rcpp::List &u_lm,
                Rcpp::List &u_glm,
                Eigen::MatrixXd &design_mat_lm,
                Eigen::MatrixXd &design_mat_glm) {
  
  int N = beta_lm.rows();
  int B = u_lm.size();
  
  Eigen::MatrixXd u_mat_lm(N, B);
  Eigen::MatrixXd u_mat_glm(N, B);
  
  for (int b = 0; b < B; ++b) {
    Eigen::VectorXd _u_lm = u_lm[b];
    Eigen::VectorXd _u_glm = u_glm[b];
    u_mat_lm.col(b) = _u_lm;
    u_mat_glm.col(b) = _u_glm;
  }
  
  Eigen::MatrixXd pred_lm = (design_mat_lm * beta_lm.transpose()) + u_mat_lm;
  Eigen::MatrixXd pred_glm = (design_mat_glm * beta_glm.transpose()) + u_mat_glm;
  pred_glm = pred_glm.unaryExpr([](double x) {return sigmoid(x); });
  
  Eigen::MatrixXd unit_preds = pred_lm.cwiseProduct(pred_glm);
  
  return Rcpp::wrap(unit_preds);
  
}

