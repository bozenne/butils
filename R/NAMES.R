### NAMES.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 23 2018 (13:09) 
## Version: 
## Last-Updated: jan  7 2020 (14:04) 
##           By: Brice Ozenne
##     Update #: 10
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

##' @title Names of an Object
##' @description Return
##' \code{colnames} for matrices,
##' \code{isS4} for S4 objects,
##' and return \code{names} otherwise.
##'
##' @param x object.
##'
##' @export
NAMES <- function(x){
    if("function" %in% class(x) || "standardGeneric" %in% class(x) ){
        return(args(x))
    }else if(isS4(x)){
        return(names(attributes(x)))
    }else if(is.matrix(x)){
        return(colnames(x))
    }else if(is.array(x)){
        return(dimnames(x))
    }else{
        return(names(x))
    }
}

##----------------------------------------------------------------------
### NAMES.R ends here
