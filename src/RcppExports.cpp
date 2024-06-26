// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// stat_resid
arma::colvec stat_resid(const arma::colvec& x, Rcpp::Nullable<NumericVector> a__, unsigned int order, const double factor, const double epsilon);
RcppExport SEXP _RcppBASTA_stat_resid(SEXP xSEXP, SEXP a__SEXP, SEXP orderSEXP, SEXP factorSEXP, SEXP epsilonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<NumericVector> >::type a__(a__SEXP);
    Rcpp::traits::input_parameter< unsigned int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< const double >::type factor(factorSEXP);
    Rcpp::traits::input_parameter< const double >::type epsilon(epsilonSEXP);
    rcpp_result_gen = Rcpp::wrap(stat_resid(x, a__, order, factor, epsilon));
    return rcpp_result_gen;
END_RCPP
}
// bin_segm
Rcpp::List bin_segm(Rcpp::List buh, double th);
RcppExport SEXP _RcppBASTA_bin_segm(SEXP buhSEXP, SEXP thSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type buh(buhSEXP);
    Rcpp::traits::input_parameter< double >::type th(thSEXP);
    rcpp_result_gen = Rcpp::wrap(bin_segm(buh, th));
    return rcpp_result_gen;
END_RCPP
}
// inner_prod_iter
arma::colvec inner_prod_iter(const arma::colvec& x);
RcppExport SEXP _RcppBASTA_inner_prod_iter(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(inner_prod_iter(x));
    return rcpp_result_gen;
END_RCPP
}
// med
double med(const arma::colvec& x);
RcppExport SEXP _RcppBASTA_med(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(med(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RcppBASTA_stat_resid", (DL_FUNC) &_RcppBASTA_stat_resid, 5},
    {"_RcppBASTA_bin_segm", (DL_FUNC) &_RcppBASTA_bin_segm, 2},
    {"_RcppBASTA_inner_prod_iter", (DL_FUNC) &_RcppBASTA_inner_prod_iter, 1},
    {"_RcppBASTA_med", (DL_FUNC) &_RcppBASTA_med, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_RcppBASTA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
