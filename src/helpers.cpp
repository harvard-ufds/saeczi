// [[Rcpp::depends(RcppEigen)]]
#include <vector>
#include <cmath>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <cli/progress.h>

double sigmoid(double x) {
  return 1.0 / (1.0 + std::exp(-x));
}

//[[Rcpp::export]]
SEXP dom_preds_calc(const Eigen::MatrixXd &beta_lm,
                    const Eigen::MatrixXd &beta_glm,
                    const int J,
                    const Eigen::MatrixXd &u_lm,
                    const Eigen::MatrixXd &u_glm,
                    const Rcpp::List &design_mats) {
  
  int B = u_lm.rows();
  
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(J, B);
  
  SEXP bar = PROTECT(cli_progress_bar(J, NULL));
  cli_progress_set_format(bar, "   - Estimating MSE {cli::pb_percent}");
  for (int j = 0; j < J; ++j) {

    Rcpp::List dmats_j = design_mats[j];
    Eigen::MatrixXd dmat_lm_j = dmats_j[0];
    Eigen::MatrixXd dmat_glm_j = dmats_j[1];
    // these are N_j x B
    Eigen::MatrixXd pred_lm_j = (dmat_lm_j * beta_lm.transpose());
    Eigen::MatrixXd pred_glm_j = (dmat_glm_j * beta_glm.transpose());
    // now need to add u's
    pred_lm_j.rowwise() += u_lm.col(j).transpose();
    pred_glm_j.rowwise() += u_glm.col(j).transpose();
    
    pred_glm_j = pred_glm_j.unaryExpr(&sigmoid);
    Eigen::MatrixXd unit_preds_j = pred_lm_j.cwiseProduct(pred_glm_j);
    
    Eigen::MatrixXd dom_preds_j = unit_preds_j.colwise().mean();
    result.block(j, 0, 1, B) = dom_preds_j;
    
    if (CLI_SHOULD_TICK) cli_progress_set(bar, j);
    
  }
  cli_progress_set_clear(bar, 1);
  cli_progress_done(bar);
  UNPROTECT(1);
  return Rcpp::wrap(result);
  
}
