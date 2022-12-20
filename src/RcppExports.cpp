// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rList
Rcpp::List rList(Rcpp::NumericMatrix m, int niter);
RcppExport SEXP _bamm_rList(SEXP mSEXP, SEXP niterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    rcpp_result_gen = Rcpp::wrap(rList(m, niter));
    return rcpp_result_gen;
END_RCPP
}
// permute_matrix
Rcpp::NumericMatrix permute_matrix(Rcpp::NumericMatrix m, int niter);
RcppExport SEXP _bamm_permute_matrix(SEXP mSEXP, SEXP niterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    rcpp_result_gen = Rcpp::wrap(permute_matrix(m, niter));
    return rcpp_result_gen;
END_RCPP
}
// null_dispersion_field_cat
Rcpp::NumericVector null_dispersion_field_cat(Rcpp::NumericMatrix dfield, Rcpp::NumericMatrix dfield_rand, double lower_interval, double upper_interval);
RcppExport SEXP _bamm_null_dispersion_field_cat(SEXP dfieldSEXP, SEXP dfield_randSEXP, SEXP lower_intervalSEXP, SEXP upper_intervalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type dfield(dfieldSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type dfield_rand(dfield_randSEXP);
    Rcpp::traits::input_parameter< double >::type lower_interval(lower_intervalSEXP);
    Rcpp::traits::input_parameter< double >::type upper_interval(upper_intervalSEXP);
    rcpp_result_gen = Rcpp::wrap(null_dispersion_field_cat(dfield, dfield_rand, lower_interval, upper_interval));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bamm_rList", (DL_FUNC) &_bamm_rList, 2},
    {"_bamm_permute_matrix", (DL_FUNC) &_bamm_permute_matrix, 2},
    {"_bamm_null_dispersion_field_cat", (DL_FUNC) &_bamm_null_dispersion_field_cat, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_bamm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}