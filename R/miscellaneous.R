### miscellaneous.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 27 2019 (09:02) 
## Version: 
## Last-Updated: sep 26 2020 (18:02) 
##           By: Brice Ozenne
##     Update #: 20
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation matrixTriangle (gdata package)
##' @title Extract or replace the upper/lower triangular portion of a matrix.
##' @description Extract or replace the upper/lower triangular portion of a matrix
##' @name matrixTriangle
##'
##' @param x [matrix]
##' @param diag [logical]
##' If \code{TRUE}, include the matrix diagonal.
##' @param byrow [logical]
##' If \code{FALSE}, return/replace elements in column-wise order.
##' If \code{TRUE}, return/replace elements in row-wise order.
##' @param value Either a single value or a vector of length equal to that of
##' the current upper/lower triangular.  Should be of a mode
##' which can be coerced to that of 'x'.
##'
##' @details Essentially a copy from the R package gdata,
##' to limit the number of dependency in butils.
##'
##' @examples
##' ## see the documentation of gdata::upperTriangle

## ** upperTriangle (code)
##' @rdname matrixTriangle
##' @export
upperTriangle <- function(x, diag=FALSE, byrow=FALSE){
    if(byrow){
        return(t(x)[rev(upper.tri(x, diag=diag))])
    } else {
        return(x[upper.tri(x, diag=diag)])
    }
}

## ** upperTriangle<- (code)
##' @rdname matrixTriangle
##' @export
`upperTriangle<-` <- function(x, diag=FALSE, byrow=FALSE, value){
    if(byrow) {
        ret <- t(x)
        ret[rev(upper.tri(x, diag=diag))] <- value
        return(t(ret))
    } else {        
        x[upper.tri(x, diag=diag)] <- value
        return(x)
    }
}

## ** lowerTriangle (code)
##' @rdname matrixTriangle
##' @export
lowerTriangle <- function(x, diag=FALSE, byrow=FALSE){
    if(byrow){
        return(t(x)[rev(lower.tri(x, diag=diag))])
    } else {
        return(x[lower.tri(x, diag=diag)])
    }
}

## ** lowerTriangle<- (code)
##' @rdname matrixTriangle
##' @export
`lowerTriangle<-` <- function(x, diag=FALSE, byrow=FALSE, value){
  if(byrow) {
    ret <- t(x)
    ret[rev(lower.tri(x, diag=diag))] <- value
    return(t(ret))
  } else {        
    x[lower.tri(x, diag=diag)] <- value
    return(x)
  }
}

## * Documentation student distribution  (LaplacesDemon package)
## copy from the LaplacesDemon package

## ** pstudent (code)
## copy from the LaplacesDemon  package
pstudent <- function (q, mu = 0, sigma = 1, nu = 10, lower.tail = TRUE, log.p = FALSE) {
    q <- as.vector(q)
    mu <- as.vector(mu)
    sigma <- as.vector(sigma)
    nu <- as.vector(nu)
    if (any(sigma <= 0)) 
        stop("The sigma parameter must be positive.")
    if (any(nu <= 0)) 
        stop("The nu parameter must be positive.")
    NN <- max(length(q), length(mu), length(sigma), length(nu))
    q <- rep(q, len = NN)
    mu <- rep(mu, len = NN)
    sigma <- rep(sigma, len = NN)
    nu <- rep(nu, len = NN)
    p <- stats::pt({q - mu}/sigma,
                   df = nu, lower.tail = lower.tail, log.p = log.p)
    temp <- which(nu > 1e+06)
    p[temp] <- stats::pnorm(q[temp], mu[temp], sigma[temp], lower.tail = lower.tail, 
                            log.p = log.p)
    return(p)
}

## ** qstudent (code)
## copy from the LaplacesDemon  package
qstudent <- function (p, mu = 0, sigma = 1, nu = 10, lower.tail = TRUE, log.p = FALSE){
    p <- as.vector(p)
    mu <- as.vector(mu)
    sigma <- as.vector(sigma)
    nu <- as.vector(nu)
    if (any(p < 0) || any(p > 1)) 
        stop("p must be in [0,1].")
    if (any(sigma <= 0)) 
        stop("The sigma parameter must be positive.")
    if (any(nu <= 0)) 
        stop("The nu parameter must be positive.")
    NN <- max(length(p), length(mu), length(sigma), length(nu))
    p <- rep(p, len = NN)
    mu <- rep(mu, len = NN)
    sigma <- rep(sigma, len = NN)
    nu <- rep(nu, len = NN)
    q <- mu + sigma * stats::qt(p,
                                df = nu,
                                lower.tail = lower.tail)
    temp <- which(nu > 1e+06)
    q[temp] <- stats::qnorm(p[temp],
                            mu[temp],
                            sigma[temp],
                            lower.tail = lower.tail, 
                            log.p = log.p)
    return(q)
}

## ** dstudent (code)
## copy from the LaplacesDemon  package
dstudent <- function(x, mu = 0, sigma = 1, nu = 10, log = FALSE){
    x <- as.vector(x)
    mu <- as.vector(mu)
    sigma <- as.vector(sigma)
    nu <- as.vector(nu)
    if (any(sigma <= 0)) 
        stop("The sigma parameter must be positive.")
    else if (any(nu <= 0)) 
        stop("The nu parameter must be positive.")
    NN <- max(length(x), length(mu), length(sigma), length(nu))
    x <- rep(x, len = NN)
    mu <- rep(mu, len = NN)
    sigma <- rep(sigma, len = NN)
    nu <- rep(nu, len = NN)
    const <- lgamma((nu + 1)/2) - lgamma(nu/2) - log(sqrt(pi * 
                                                          nu) * sigma)
    dens <- const + log((1 + (1/nu) * ((x - mu)/sigma)^2)^(-(nu + 
                                                             1)/2))
    if (log == FALSE) 
        dens <- exp(dens)
    return(dens)
}


######################################################################
### miscellaneous.R ends here



