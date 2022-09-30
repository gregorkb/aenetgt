// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// EYgibbs
Rcpp::NumericVector EYgibbs(int N, NumericVector p, NumericMatrix Y, NumericMatrix Z, NumericVector se, NumericVector sp, int na, int GI);
RcppExport SEXP _aenetgt_EYgibbs(SEXP NSEXP, SEXP pSEXP, SEXP YSEXP, SEXP ZSEXP, SEXP seSEXP, SEXP spSEXP, SEXP naSEXP, SEXP GISEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type se(seSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sp(spSEXP);
    Rcpp::traits::input_parameter< int >::type na(naSEXP);
    Rcpp::traits::input_parameter< int >::type GI(GISEXP);
    rcpp_result_gen = Rcpp::wrap(EYgibbs(N, p, Y, Z, se, sp, na, GI));
    return rcpp_result_gen;
END_RCPP
}
// EYiYjgibbs_slow
Rcpp::NumericMatrix EYiYjgibbs_slow(int N, NumericVector p, NumericMatrix Y, NumericMatrix Z, NumericVector se, NumericVector sp, int na, int GI);
RcppExport SEXP _aenetgt_EYiYjgibbs_slow(SEXP NSEXP, SEXP pSEXP, SEXP YSEXP, SEXP ZSEXP, SEXP seSEXP, SEXP spSEXP, SEXP naSEXP, SEXP GISEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type se(seSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sp(spSEXP);
    Rcpp::traits::input_parameter< int >::type na(naSEXP);
    Rcpp::traits::input_parameter< int >::type GI(GISEXP);
    rcpp_result_gen = Rcpp::wrap(EYiYjgibbs_slow(N, p, Y, Z, se, sp, na, GI));
    return rcpp_result_gen;
END_RCPP
}
// CovYiYjgibbs
Rcpp::NumericMatrix CovYiYjgibbs(int N, NumericVector p, NumericMatrix Y, NumericMatrix Z, NumericMatrix W, NumericVector se, NumericVector sp, NumericVector EY, int na, int GI);
RcppExport SEXP _aenetgt_CovYiYjgibbs(SEXP NSEXP, SEXP pSEXP, SEXP YSEXP, SEXP ZSEXP, SEXP WSEXP, SEXP seSEXP, SEXP spSEXP, SEXP EYSEXP, SEXP naSEXP, SEXP GISEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type W(WSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type se(seSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sp(spSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type EY(EYSEXP);
    Rcpp::traits::input_parameter< int >::type na(naSEXP);
    Rcpp::traits::input_parameter< int >::type GI(GISEXP);
    rcpp_result_gen = Rcpp::wrap(CovYiYjgibbs(N, p, Y, Z, W, se, sp, EY, na, GI));
    return rcpp_result_gen;
END_RCPP
}
// logistic_enet
Rcpp::List logistic_enet(Rcpp::NumericVector Yr, Rcpp::NumericMatrix Xr, float lambda, Rcpp::NumericVector gammar, float theta, float delta);
RcppExport SEXP _aenetgt_logistic_enet(SEXP YrSEXP, SEXP XrSEXP, SEXP lambdaSEXP, SEXP gammarSEXP, SEXP thetaSEXP, SEXP deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Yr(YrSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Xr(XrSEXP);
    Rcpp::traits::input_parameter< float >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type gammar(gammarSEXP);
    Rcpp::traits::input_parameter< float >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< float >::type delta(deltaSEXP);
    rcpp_result_gen = Rcpp::wrap(logistic_enet(Yr, Xr, lambda, gammar, theta, delta));
    return rcpp_result_gen;
END_RCPP
}
// llj_array
Rcpp::List llj_array(Rcpp::IntegerVector Zjr, Rcpp::IntegerVector Zjc, Rcpp::IntegerVector Yji, Rcpp::IntegerVector whichjretest, Rcpp::NumericVector pxji, Rcpp::NumericVector Se, Rcpp::NumericVector Sp, int B);
RcppExport SEXP _aenetgt_llj_array(SEXP ZjrSEXP, SEXP ZjcSEXP, SEXP YjiSEXP, SEXP whichjretestSEXP, SEXP pxjiSEXP, SEXP SeSEXP, SEXP SpSEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Zjr(ZjrSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Zjc(ZjcSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Yji(YjiSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type whichjretest(whichjretestSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pxji(pxjiSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Se(SeSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Sp(SpSEXP);
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(llj_array(Zjr, Zjc, Yji, whichjretest, pxji, Se, Sp, B));
    return rcpp_result_gen;
END_RCPP
}
// all_binary_sequences
arma::mat all_binary_sequences(int a);
RcppExport SEXP _aenetgt_all_binary_sequences(SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(all_binary_sequences(a));
    return rcpp_result_gen;
END_RCPP
}
// EYexact
arma::colvec EYexact(IntegerMatrix Z, IntegerMatrix Y, NumericMatrix X, NumericVector b, NumericVector Se, NumericVector Sp);
RcppExport SEXP _aenetgt_EYexact(SEXP ZSEXP, SEXP YSEXP, SEXP XSEXP, SEXP bSEXP, SEXP SeSEXP, SEXP SpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Se(SeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Sp(SpSEXP);
    rcpp_result_gen = Rcpp::wrap(EYexact(Z, Y, X, b, Se, Sp));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_aenetgt_EYgibbs", (DL_FUNC) &_aenetgt_EYgibbs, 8},
    {"_aenetgt_EYiYjgibbs_slow", (DL_FUNC) &_aenetgt_EYiYjgibbs_slow, 8},
    {"_aenetgt_CovYiYjgibbs", (DL_FUNC) &_aenetgt_CovYiYjgibbs, 10},
    {"_aenetgt_logistic_enet", (DL_FUNC) &_aenetgt_logistic_enet, 6},
    {"_aenetgt_llj_array", (DL_FUNC) &_aenetgt_llj_array, 8},
    {"_aenetgt_all_binary_sequences", (DL_FUNC) &_aenetgt_all_binary_sequences, 1},
    {"_aenetgt_EYexact", (DL_FUNC) &_aenetgt_EYexact, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_aenetgt(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
