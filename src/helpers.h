#ifndef helpers_H
#define helpers_H

#include <Rcpp.h>
#include <RcppEigen.h>

double sigmoid(double x);

void preds_calc(Eigen::MatrixXd& result,
                const Eigen::MatrixXd& beta_lm,
                const Eigen::MatrixXd& beta_glm,
                const Eigen::MatrixXd& dmat_lm,
                const Eigen::MatrixXd& dmat_glm,
                const Eigen::MatrixXd& u_lm,
                const Eigen::MatrixXd& u_glm,
                int j,
                int B,
                std::string estimand);

SEXP generate_preds(const Eigen::MatrixXd& beta_lm,
                    const Eigen::MatrixXd& beta_glm,
                    const Eigen::MatrixXd& u_lm,
                    const Eigen::MatrixXd& u_glm,
                    const Rcpp::List& design_mats,
                    int J,
                    std::string estimand);

#endif // helpers_H