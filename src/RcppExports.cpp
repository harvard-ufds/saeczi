// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// generate_preds
SEXP generate_preds(const Eigen::MatrixXd& beta_lm, const Eigen::MatrixXd& beta_glm, const Eigen::MatrixXd& u_lm, const Eigen::MatrixXd& u_glm, const Rcpp::List& design_mats, int J, std::string estimand, Rcpp::Nullable<Rcpp::Function> inv);
RcppExport SEXP _saeczi_generate_preds(SEXP beta_lmSEXP, SEXP beta_glmSEXP, SEXP u_lmSEXP, SEXP u_glmSEXP, SEXP design_matsSEXP, SEXP JSEXP, SEXP estimandSEXP, SEXP invSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type beta_lm(beta_lmSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type beta_glm(beta_glmSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type u_lm(u_lmSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type u_glm(u_glmSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type design_mats(design_matsSEXP);
    Rcpp::traits::input_parameter< int >::type J(JSEXP);
    Rcpp::traits::input_parameter< std::string >::type estimand(estimandSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::Function> >::type inv(invSEXP);
    rcpp_result_gen = Rcpp::wrap(generate_preds(beta_lm, beta_glm, u_lm, u_glm, design_mats, J, estimand, inv));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_saeczi_generate_preds", (DL_FUNC) &_saeczi_generate_preds, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_saeczi(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
