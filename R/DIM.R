### DIM.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 22 2018 (15:06) 
## Version: 
## Last-Updated: jun 22 2018 (15:09) 
##           By: Brice Ozenne
##     Update #: 1
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

##' @title Dimensions of an Object
##' @description Return \code{length} for vector and list and use \code{dim} otherwise.
##'
##' @param x object.
##'
##' @export
DIM <- function(x){
    if(is.list(x)||is.vector(x)){
        return(length(x))
    }else{
        retrun(dim(x))
    }
}

######################################################################
### DIM.R ends here
