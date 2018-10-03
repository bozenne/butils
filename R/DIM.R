### DIM.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 22 2018 (15:06) 
## Version: 
## Last-Updated: okt  3 2018 (14:29) 
##           By: Brice Ozenne
##     Update #: 7
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
    } else if(data.table::is.data.table(x)){
        return(dim(x))
    } else if(is.data.frame(x)){
        return(dim(x))
    } else if(is.list(x)||is.vector(x)){
        return(length(x))
    } else {
        return(dim(x))
    }
}

######################################################################
### DIM.R ends here
