### WaldTest.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 23 2017 (14:57) 
## Version: 
## Last-Updated: nov 24 2017 (09:50) 
##           By: Brice Ozenne
##     Update #: 68
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * WaldTest - documentation
#' @title Testing independent linear hypotheses using a Wald test
#' @description Use a contrast matrix to test several linear hypotheses on the coefficients using a Wald test.
#' The hypotheses are assumed to be independent of each other.
#' @name WaldTest
#' 
#' @param object a model. Currently only support \code{lme} objects.
#' @param C a contrast matrix.
#' Number of rows is the number of hypotheses.
#' Number of columns should match the number of coefficients in the model.
#' @param b a vector such that the hypothesis to test is: \eqn{C * \beta = b} where \eqn{\beta} are the model coefficients.
#' @param df [optional] the degree of freedom associated to the variance of each coefficient.
#' @param ... not used.
#' 
#' @details
#' Denoting \eqn{\beta} the estimated model coefficients,
#' \eqn{Sigma} their estimated variance covariance matrix,
#' and \eqn{t()} the transpose operator, this function computes:
#' \deqn{Cb = C * \beta}
#' \deqn{CSC = C * \Sigma * t(C)}
#' \deqn{t = \sqrt{Cb/diag(CSC)}}
#'
#' and compute the p.value using the Gaussian distribution (\code{df=NULL}) or a student's t-distribution.
#' In such a case the degrees of freedom are computed using:
#' \deqn{Cdf = C * df}
#' This formula is probably a very crude approximation to the appropriate degrees of freedom of \eqn{diage(CSC)} (assuming they exists).
#'
#' @export
`WaldTest` <-
    function(object, ...) UseMethod("WaldTest")

## * WaldTest.lme
#' @rdname WaldTest
#' @export
WaldTest.lme <- function(object, C, b = rep(0,NROW(C)), df = NULL){

    beta <- fixef(object)
    Sigma <- stats::vcov(object)
    p <- length(beta)
    n <- object$dims$N
    Sigma.corrected <- Sigma * n/(n-p)
    if(is.null(df)){
        df <- rep(n-p,p)
    }
    
    .WaldTest(beta = beta, Sigma = Sigma.corrected, C = C, b = b, df = df)
}

## * .WaldTest
.WaldTest <- function(beta, Sigma, C, b, df){
### ** prepare
    name.coef <- names(beta)
    if(is.null(name.coef)){
        stop("the coefficients must be named")
    }
    n.coef <- length(name.coef)
    if(NCOL(C) != n.coef){
        stop("the number of columns of the constrast matrix should be ",n.coef,"\n")
    }
    if(NROW(C) != length(b)){
        stop("the number of rows of the constrast matrix should match the length of argument \'b\' \n")
    }
    
### ** compute statistic and p.value
    Cbeta <- as.vector(C %*% beta)    
    CSigmaC <- C %*% Sigma %*% t(C)
    ## average of df over parameters
    if(!is.null(df)){
        Cdf <- as.vector(C %*% df)/rowSums(C)
    }else{
        Cdf <- rep(NA, length(Cbeta) )
    }

    diag.CSigmaC <- sqrt(diag(C %*% Sigma %*% t(C)))
    t.stat <- Cbeta/diag.CSigmaC
    if(is.null(df)){
        p.value <- 2*(1-stats::pnorm(abs(t.stat)))
    }else{
        p.value <- 2*(1-stats::pt(abs(t.stat), df = Cdf))
    }
    signif <- sapply(p.value, function(iP){
        if(iP<0.001){
            return("***")
        }else if(iP<0.01){
            return("**")
        }else if(iP<0.05){
            return("*")
        }else if(iP<0.1){
            return(".")
        }else{
            return("")
        }
    })
    out <- data.frame("Estimate" = Cbeta,
                      "Std.Error" = diag.CSigmaC,
                      "t.value" =  t.stat,
                      "df" = Cdf,
                      "p.value" = p.value,
                      "significance" = signif)

    return(out)
}


##----------------------------------------------------------------------
### WaldTest.R ends here
