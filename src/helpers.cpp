// [[Rcpp::depends(RcppEigen)]]
#include <vector>
#include <unordered_map>
#include <cmath>
#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Eigen;
using namespace Rcpp;

Eigen::VectorXd match_val(const Rcpp::CharacterVector &names,
                          const Rcpp::NumericVector &values,
                          const Rcpp::CharacterVector &input) {
  
  std::unordered_map<Rcpp::String, double> nameMap;
  
  for (size_t i = 0; i < names.size(); ++i) {
    nameMap[names[i]] = values[i];
  }
  
  std::vector<double> result;
  result.reserve(input.size());
  
  for (const auto& i : input) {
    
    auto it = nameMap.find(i);
    
    if (it != nameMap.end()) {
      result.push_back(it->second);
    } else {
      // this is like allow.new.levels in R
      result.push_back(0);
    }
  }
  
  Eigen::Map<Eigen::VectorXd> res_eigen(result.data(), result.size());
  return res_eigen;
}

void add_u(const Rcpp::CharacterVector &names,
           const Rcpp::List &u,
           const Rcpp::CharacterVector &input,
           Eigen::MatrixXd &lin_comp_lm,
           Eigen::MatrixXd &lin_comp_glm,
           const int B) {
  
  for(int b = 0; b < B; ++b) {
    
    Rcpp::List u_b = u[b];
    Rcpp::NumericVector u_b_lm = u_b[0];
    Rcpp::NumericVector u_b_glm = u_b[1];
    Eigen::VectorXd u_b_lm_match = match_val(names, u_b_lm, input);
    Eigen::VectorXd u_b_glm_match = match_val(names, u_b_glm, input);
    
    lin_comp_lm.col(b) = lin_comp_lm.col(b) + u_b_lm_match;
    lin_comp_glm.col(b) = lin_comp_glm.col(b) + u_b_glm_match;
    
  }
  
}

double sigmoid(double x) {
  return 1.0 / (1.0 + std::exp(-x));
}


//[[Rcpp::export]]
SEXP predict_cpp(const Eigen::MatrixXd &beta_lm,
                 const Eigen::MatrixXd &beta_glm,
                 const Rcpp::CharacterVector &names,
                 const Rcpp::CharacterVector &dom_input,
                 const Rcpp::List &u,
                 const Eigen::MatrixXd &design_mat_lm,
                 const Eigen::MatrixXd &design_mat_glm) {

  int B = u.size();
  
  Eigen::MatrixXd pred_lm = (design_mat_lm * beta_lm.transpose());
  Eigen::MatrixXd pred_glm = (design_mat_glm * beta_glm.transpose());
  
  add_u(names, u, dom_input, pred_lm, pred_glm, B);
  
  pred_glm = pred_glm.array().unaryExpr([](double x) {return sigmoid(x); });
  
  Eigen::MatrixXd unit_preds = pred_lm.cwiseProduct(pred_glm);
  
  return Rcpp::wrap(unit_preds);
  
}

