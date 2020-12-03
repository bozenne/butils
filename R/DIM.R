### DIM.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 22 2018 (15:06) 
## Version: 
## Last-Updated: nov 26 2020 (18:27) 
##           By: Brice Ozenne
##     Update #: 8
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
    } else if(is.data.frame(x)||is.matrix(x)){
        return(dim(x))
    } else if(is.list(x)||is.vector(x)||is.null(dim(x))){
        return(length(x))
    } else {
        return(dim(x))
    }
}

######################################################################
### DIM.R ends here
