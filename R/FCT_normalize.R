### normalize.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: apr 25 2017 (14:13) 
## Version: 
## last-updated: apr 25 2017 (15:13) 
##           By: Brice Ozenne
##     Update #: 16
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

EDF <- function(X,x,n){
    mean(X<=x)-1/(2*n)
}
erfinv <- function (x){ # from http://stackoverflow.com/questions/29067916/r-error-function-erfz
    qnorm((1 + x)/2)/sqrt(2)
}

#' @title Transform a variable to obtain approximate normality
#' @description (Slow) implementation of the transformation described in Albada, 2007.
#' 
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
#' shapiro.test(Xnorm)
#' # plot(X, Xnorm)
#'
#' X.NA <- c(NA,X,NA)
#' Xnorm.NA <- normalize(X.NA, na.rm = TRUE)
#' 
#' ## gamma distribution
#' X <- rgamma(n, shape = 1)
#' shapiro.test(X)
#' # hist(X)
#' 
#' Xnorm <- normalize(X)
#' shapiro.test(Xnorm)
#' # hist(Xnorm)
#' # plot(X,Xnorm)
#' @export
normalize <- function(X, na.rm = FALSE){
  
  if(!is.numeric(X)){
    stop("\'X\' must be a numeric vector \n")
  }
  
  if(any(is.na(X))){
    if(na.rm){
      Xsave <- X
      X <- as.numeric(na.omit(X))
      test.NA <- TRUE
    }else{
      stop("\'X\' contains NA \n",
           "set na.rm=TRUE to ignore them \n")
    }
  }else{
    test.NA <- FALSE
  }
  
    n <- length(X)
    X <- sapply(X, function(x){
        sqrt(2)*erfinv(2*EDF(X,x,n)-1)
    })
    
    if(test.NA){
      Xsave[!is.na(Xsave)] <- X
      X <- Xsave
    }
    
    return(X)
}



#----------------------------------------------------------------------
### normalize.R ends here
