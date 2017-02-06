// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/butils.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

// colCumSum
NumericMatrix colCumSum(NumericMatrix x);
RcppExport SEXP butils_colCumSum(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(colCumSum(x));
    return rcpp_result_gen;
END_RCPP
}
// hello0
void hello0();
static SEXP butils_hello0_try() {
BEGIN_RCPP
    hello0();
    return R_NilValue;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP butils_hello0() {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(butils_hello0_try());
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// hello1
arma::mat hello1(arma::mat A);
static SEXP butils_hello1_try(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(hello1(A));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP butils_hello1(SEXP ASEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(butils_hello1_try(ASEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// hello2
arma::mat hello2(arma::mat A, SEXP B);
static SEXP butils_hello2_try(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< SEXP >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(hello2(A, B));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP butils_hello2(SEXP ASEXP, SEXP BSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(butils_hello2_try(ASEXP, BSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}

// validate (ensure exported C++ functions exist before calling them)
static int butils_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("void(*hello0)()");
        signatures.insert("arma::mat(*hello1)(arma::mat)");
        signatures.insert("arma::mat(*hello2)(arma::mat,SEXP)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP butils_RcppExport_registerCCallable() { 
    R_RegisterCCallable("butils", "butils_hello0", (DL_FUNC)butils_hello0_try);
    R_RegisterCCallable("butils", "butils_hello1", (DL_FUNC)butils_hello1_try);
    R_RegisterCCallable("butils", "butils_hello2", (DL_FUNC)butils_hello2_try);
    R_RegisterCCallable("butils", "butils_RcppExport_validate", (DL_FUNC)butils_RcppExport_validate);
    return R_NilValue;
}
