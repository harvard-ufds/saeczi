// [[Rcpp::depends(RcppEigen)]]
#include <vector>
#include <unordered_map>
#include <cmath>
#include <Rcpp.h>
#include <RcppEigen.h>

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

void squish_rows(Eigen::MatrixXd &unit_preds_mat,
                 Eigen::MatrixXd &result_mat,
                 std::vector<std::string> &unique_doms,
                 std::unordered_map<std::string, std::vector<int>> &doms_map) {
  
  int n_rows = unique_doms.size();
  int n_cols = unit_preds_mat.cols();
  
  for (int i = 0; i < n_rows; ++i) {
    
    const std::vector<int>& dom_ids = doms_map[unique_doms[i]];
    
    // subset the matrix to only include rows from dom_ids
    Eigen::MatrixXd subset(dom_ids.size(), n_cols);
    for (int d = 0; d < dom_ids.size(); ++d) {
      subset.row(d) = unit_preds_mat.row(dom_ids[d]);
    }

    Eigen::MatrixXd squished_subset = subset.colwise().mean();
    
    result_mat.block(i, 0, 1, n_cols) = squished_subset;
    
  }
  
}


//[[Rcpp::export]]
SEXP dom_preds_calc(const Eigen::MatrixXd &beta_lm,
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
  
  std::vector<std::string> dom_input_std(dom_input.size());
  for (int i = 0; i < dom_input.size(); ++i) {
    dom_input_std[i] = Rcpp::as<std::string>(dom_input[i]);
  }

  std::sort(dom_input_std.begin(), dom_input_std.end());
  auto unique_doms = std::unique(dom_input_std.begin(), dom_input_std.end());
  dom_input_std.erase(unique_doms, dom_input_std.end());

  int n_doms = dom_input_std.size();

  Eigen::MatrixXd result(n_doms, unit_preds.cols());

  std::unordered_map<std::string, std::vector<int>> dom_id_map;
  for (int i = 0; i < dom_input.size(); ++i) {
    dom_id_map[dom_input_std[i]].push_back(i);
  }

  squish_rows(unit_preds,
              result,
              dom_input_std,
              dom_id_map);

  Rcpp::CharacterVector row_doms = Rcpp::wrap(dom_input_std);
  Rcpp::List out = Rcpp::List::create(result, row_doms);
  
  return Rcpp::wrap(out);
  
}

