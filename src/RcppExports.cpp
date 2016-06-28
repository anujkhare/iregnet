// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// compute_grad_response_cpp
Rcpp::List compute_grad_response_cpp(Rcpp::NumericVector y_l, Rcpp::NumericVector y_r, Rcpp::NumericVector eta, double scale, Rcpp::IntegerVector censoring_type, Rcpp::String family);
RcppExport SEXP iregnet_compute_grad_response_cpp(SEXP y_lSEXP, SEXP y_rSEXP, SEXP etaSEXP, SEXP scaleSEXP, SEXP censoring_typeSEXP, SEXP familySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y_l(y_lSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y_r(y_rSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type censoring_type(censoring_typeSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type family(familySEXP);
    __result = Rcpp::wrap(compute_grad_response_cpp(y_l, y_r, eta, scale, censoring_type, family));
    return __result;
END_RCPP
}
// compute_densities
Rcpp::NumericVector compute_densities(Rcpp::NumericVector z, int j, Rcpp::String family);
RcppExport SEXP iregnet_compute_densities(SEXP zSEXP, SEXP jSEXP, SEXP familySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type family(familySEXP);
    __result = Rcpp::wrap(compute_densities(z, j, family));
    return __result;
END_RCPP
}
// fit_cpp
Rcpp::List fit_cpp(Rcpp::NumericMatrix X, Rcpp::NumericMatrix y, Rcpp::String family, double alpha, bool intercept, double scale, bool standardize, double max_iter, double threshold, int num_lambda, double eps_lambda, int flag_debug);
RcppExport SEXP iregnet_fit_cpp(SEXP XSEXP, SEXP ySEXP, SEXP familySEXP, SEXP alphaSEXP, SEXP interceptSEXP, SEXP scaleSEXP, SEXP standardizeSEXP, SEXP max_iterSEXP, SEXP thresholdSEXP, SEXP num_lambdaSEXP, SEXP eps_lambdaSEXP, SEXP flag_debugSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type family(familySEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< bool >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< bool >::type standardize(standardizeSEXP);
    Rcpp::traits::input_parameter< double >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< int >::type num_lambda(num_lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type eps_lambda(eps_lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type flag_debug(flag_debugSEXP);
    __result = Rcpp::wrap(fit_cpp(X, y, family, alpha, intercept, scale, standardize, max_iter, threshold, num_lambda, eps_lambda, flag_debug));
    return __result;
END_RCPP
}
