// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// unit_preds_calc
SEXP unit_preds_calc(const Eigen::MatrixXd& beta_lm, const Eigen::MatrixXd& beta_glm, const Rcpp::CharacterVector& names, const Rcpp::CharacterVector& dom_input, const Rcpp::List& u, const Eigen::MatrixXd& design_mat_lm, const Eigen::MatrixXd& design_mat_glm);
RcppExport SEXP _saeczi_unit_preds_calc(SEXP beta_lmSEXP, SEXP beta_glmSEXP, SEXP namesSEXP, SEXP dom_inputSEXP, SEXP uSEXP, SEXP design_mat_lmSEXP, SEXP design_mat_glmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type beta_lm(beta_lmSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type beta_glm(beta_glmSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type names(namesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type dom_input(dom_inputSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type design_mat_lm(design_mat_lmSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type design_mat_glm(design_mat_glmSEXP);
    rcpp_result_gen = Rcpp::wrap(unit_preds_calc(beta_lm, beta_glm, names, dom_input, u, design_mat_lm, design_mat_glm));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_saeczi_unit_preds_calc", (DL_FUNC) &_saeczi_unit_preds_calc, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_saeczi(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
