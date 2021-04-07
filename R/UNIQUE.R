### UNIQUE.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Apr  7 2021 (11:09) 
## Version: 
## Last-Updated: Apr  7 2021 (11:12) 
##           By: Brice Ozenne
##     Update #: 3
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

##' @title Unique Elements of an Object
##' @description Return \code{levels(droplevels(.))} for factor variables and otherwise \code{unique(.)}
##'
##' @param x object.
##'
##' @export
UNIQUE <- function(x){
    if(is.factor(x)){
        return(levels(droplevels(x)))
    }else{
        return(unique(x))
    }
}


##----------------------------------------------------------------------
### UNIQUE.R ends here
