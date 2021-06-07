// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// colSumsSub
arma::rowvec colSumsSub(const arma::mat& mat, const arma::uvec& rowIDs);
RcppExport SEXP _MSETools_colSumsSub(SEXP matSEXP, SEXP rowIDsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type rowIDs(rowIDsSEXP);
    rcpp_result_gen = Rcpp::wrap(colSumsSub(mat, rowIDs));
    return rcpp_result_gen;
END_RCPP
}
// dec
double dec(const Rcpp::IntegerVector I, int p);
RcppExport SEXP _MSETools_dec(SEXP ISEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type I(ISEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(dec(I, p));
    return rcpp_result_gen;
END_RCPP
}
// decClique
arma::uvec decClique(const Rcpp::List cliques, int p);
RcppExport SEXP _MSETools_decClique(SEXP cliquesSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type cliques(cliquesSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(decClique(cliques, p));
    return rcpp_result_gen;
END_RCPP
}
// computeML
void computeML(arma::mat& inPlace, int j, const arma::mat& compMat, const Rcpp::List& C, const Rcpp::List& S, const arma::rowvec& denominator, int p);
RcppExport SEXP _MSETools_computeML(SEXP inPlaceSEXP, SEXP jSEXP, SEXP compMatSEXP, SEXP CSEXP, SEXP SSEXP, SEXP denominatorSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type inPlace(inPlaceSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type compMat(compMatSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type C(CSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type denominator(denominatorSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    computeML(inPlace, j, compMat, C, S, denominator, p);
    return R_NilValue;
END_RCPP
}
// computeLogPostProbs
arma::mat computeLogPostProbs(const arma::mat& compMat, const Rcpp::List& graphs, const arma::rowvec& denominator, int p);
RcppExport SEXP _MSETools_computeLogPostProbs(SEXP compMatSEXP, SEXP graphsSEXP, SEXP denominatorSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type compMat(compMatSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type graphs(graphsSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type denominator(denominatorSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(computeLogPostProbs(compMat, graphs, denominator, p));
    return rcpp_result_gen;
END_RCPP
}
// rowAdd
void rowAdd(arma::mat& mat, const arma::rowvec& v);
RcppExport SEXP _MSETools_rowAdd(SEXP matSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type v(vSEXP);
    rowAdd(mat, v);
    return R_NilValue;
END_RCPP
}
// colAdd
void colAdd(arma::mat& mat, const arma::colvec& v);
RcppExport SEXP _MSETools_colAdd(SEXP matSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type v(vSEXP);
    colAdd(mat, v);
    return R_NilValue;
END_RCPP
}
// expNormalize
void expNormalize(arma::mat& mat);
RcppExport SEXP _MSETools_expNormalize(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type mat(matSEXP);
    expNormalize(mat);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MSETools_colSumsSub", (DL_FUNC) &_MSETools_colSumsSub, 2},
    {"_MSETools_dec", (DL_FUNC) &_MSETools_dec, 2},
    {"_MSETools_decClique", (DL_FUNC) &_MSETools_decClique, 2},
    {"_MSETools_computeML", (DL_FUNC) &_MSETools_computeML, 7},
    {"_MSETools_computeLogPostProbs", (DL_FUNC) &_MSETools_computeLogPostProbs, 4},
    {"_MSETools_rowAdd", (DL_FUNC) &_MSETools_rowAdd, 2},
    {"_MSETools_colAdd", (DL_FUNC) &_MSETools_colAdd, 2},
    {"_MSETools_expNormalize", (DL_FUNC) &_MSETools_expNormalize, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_MSETools(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
