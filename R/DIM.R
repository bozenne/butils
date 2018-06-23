### DIM.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 22 2018 (15:06) 
## Version: 
## Last-Updated: jun 23 2018 (13:14) 
##           By: Brice Ozenne
##     Update #: 5
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

##' @title Dimensions of an Object
##' @description Return \code{length} for vectors, lists, S4 objects,
##' and use \code{dim} otherwise.
##'
##' @param x object.
##'
##' @export
DIM <- function(x){
    if(isS4(x)){
        length(attributes(x))
    }else if(is.list(x)||is.vector(x)){
        return(length(x))
    }else{
        retrun(dim(x))
    }
}

######################################################################
### DIM.R ends here
