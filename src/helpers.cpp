// [[Rcpp::depends(RcppEigen)]]
#include <vector>
#include <unordered_map>
#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Eigen;
using namespace Rcpp;

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
  
  //Eigen::MatrixXd* u_lm_ptr = &u_mat_lm;
  //Eigen::MatrixXd* u_glm_ptr = &u_mat_glm;
  
  
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

//[[Rcpp::export]]
SEXP match_val(const Rcpp::CharacterVector& names,
               const Rcpp::NumericVector& values,
               const Rcpp::CharacterVector& input) {
  
  std::unordered_map<Rcpp::String, double> nameMap;
  
  for (size_t i = 0; i < names.size(); ++i) {
    nameMap[names[i]] = values[i];
  }
  
  std::vector<double> result;
  result.reserve(input.size());
  
  for (const auto& inp : input) {

    auto it = nameMap.find(inp);
    
    if (it != nameMap.end()) {
      result.push_back(it->second);
    } else {
      // this is like allow.new.levels in R
      result.push_back(0);
    }
  }
  
  return Rcpp::wrap(result);
}
