// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_butils_RCPPEXPORTS_H_GEN_
#define RCPP_butils_RCPPEXPORTS_H_GEN_

#include <RcppArmadillo.h>
#include <Rcpp.h>

namespace butils {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("butils", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("butils", "butils_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in butils");
            }
        }
    }

    inline void hello0() {
        typedef SEXP(*Ptr_hello0)();
        static Ptr_hello0 p_hello0 = NULL;
        if (p_hello0 == NULL) {
            validateSignature("void(*hello0)()");
            p_hello0 = (Ptr_hello0)R_GetCCallable("butils", "butils_hello0");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_hello0();
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
    }

    inline arma::mat hello1(arma::mat A) {
        typedef SEXP(*Ptr_hello1)(SEXP);
        static Ptr_hello1 p_hello1 = NULL;
        if (p_hello1 == NULL) {
            validateSignature("arma::mat(*hello1)(arma::mat)");
            p_hello1 = (Ptr_hello1)R_GetCCallable("butils", "butils_hello1");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_hello1(Rcpp::wrap(A));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

    inline arma::mat hello2(arma::mat A, SEXP B) {
        typedef SEXP(*Ptr_hello2)(SEXP,SEXP);
        static Ptr_hello2 p_hello2 = NULL;
        if (p_hello2 == NULL) {
            validateSignature("arma::mat(*hello2)(arma::mat,SEXP)");
            p_hello2 = (Ptr_hello2)R_GetCCallable("butils", "butils_hello2");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_hello2(Rcpp::wrap(A), Rcpp::wrap(B));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

}

#endif // RCPP_butils_RCPPEXPORTS_H_GEN_
