### extractData.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec  8 2019 (22:28) 
## Version: 
## Last-Updated: jan 23 2020 (09:28) 
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

## * method extractData.lmer
#' @rdname extractData
#' @export
extractData.merMod <- function(object, design.matrix = FALSE, as.data.frame = TRUE){
    ## ** check arguments
    validLogical(design.matrix, valid.length = 1)
    validLogical(as.data.frame, valid.length = 1)

    ## ** extract data
    if(design.matrix){
        data <- try(model.matrix(object), silent = TRUE)
    }else{
        data <- evalInParentEnv(object@call$data)
        
        if("function" %in% class(data)){
            stop("data has the same name as a function \n",
                 "consider renaming data before generating object \n")
        }
        if(!inherits(data, "data.frame")){
            stop("Could not extract the data from the model \n")
        } 
    }

    ## ** normalize data
    if(as.data.frame){
        data <- as.data.frame(data)        
    }

    ## ** export
    return(data)
    
}

##----------------------------------------------------------------------
### extractData.R ends here
