### normalize.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: apr 25 2017 (14:13) 
## Version: 
## last-updated: apr 25 2017 (14:19) 
##           By: Brice Ozenne
##     Update #: 7
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

EDF <- function(x){
    mean(X<=x)
}

erfinv <- function (x){
    qnorm((1 + x)/2)/sqrt(2)
}

#' @title Transform a variable to obtain approximate normality
#' @param X the vector of values
#' 
#' @references Albada et al. Transformation of arbitrary distributions to the normal distribution with application to EEG test-retest reliability (2007, journal of Neuroscience Methods)
#' @examples
#' 
#' n <- 1000
#'
#' ## normal distribution ##
#' X <- rnorm(n)
#' Xnorm <- normalize(X)
#' shapiro.test(Xnorm[!is.infinite(Xnorm)])
#' # plot(X, Xnorm)
#'
#' ## gamma distribution
#' X <- rgamma(n, shape = 1)
#' shapiro.test(X)
#' # hist(X)
#' 
#' Xnorm <- normalize(X)
#' shapiro.test(Xnorm[!is.infinite(Xnorm)])
#' # hist(Xnorm[!is.infinite(Xnorm)])
#' plot(X[!is.infinite(Xnorm)],Xnorm[!is.infinite(Xnorm)])
#' @export
normalize <- function(X){
    sapply(X, function(x){
        sqrt(2)*erfinv(2*EDF(x)-1)
    })
}



#----------------------------------------------------------------------
### normalize.R ends here
