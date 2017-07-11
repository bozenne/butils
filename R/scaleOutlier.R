### scaleOutlier.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 25 2017 (18:49) 
## Version: 
## last-updated: maj 25 2017 (22:20) 
##           By: Brice Ozenne
##     Update #: 2
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#' @title Robust scaling
#' @description Robust scaling of a numeric vector 
#'
#' @param x a numeric vector
#' @param center the method used to assess the "average" value
#' @param scale the method used to assess the dispersion of the data
#' @param method a method taking the dataset as its first argument and returning first the center and then the scale value. Disregarded if arguments center and scale are not null.
#' @param na.rm should missing values be omitted
#' @param noScaleIf0 If \code{TRUE} the variable will not be scaled if its dispersion is null.
#' @param ... additional arguments passed to method
#' 
#' @return the scaled vector
#' 
#' @examples 
#' n <- 1e3
#' scaleOutlier(rnorm(n))
#' 
#' @export
scaleOutlier <- function(x, center = "median", scale = "mad", method, 
                         na.rm = FALSE, noScaleIf0 = FALSE, ...){
  
  if(na.rm == TRUE){
    xx <- na.omit(x)
  }else{
    xx <- x
  }
  
  if(is.null(center) && is.null(scale) && !missing(method)){
    
    res.method <- method(xx, ...)
    centerR <- res.method[1]
    scaleR <- res.method[2]
    
  }else{
    
    if(is.character(center)){
      centerR <- do.call(center, args = list(xx, ...))
    }else if(is.numeric(center)){
      centerR <- center
    }else if(identical(center, FALSE)){
      centerR <- FALSE
    }else{
      stop("scaleOutlier: unknown type of center parameter \n")
    }
    
    if(is.character(scale)){
      scaleR <- do.call(scale, args = list(xx, ...))
    }else if(is.numeric(scale)){
      scaleR <- scale
    }else if(identical(scale, FALSE)){
      scaleR <- FALSE
    }else{
      stop("scaleOutlier: unknown type of scale parameter \n")
    }
    
  }
  
  if(is.infinite(centerR)){ 
    stop("scaleOutlier: infinite center parameter \n")
  }
  if(is.infinite(scaleR)){ 
    stop("scaleOutlier: infinite scale parameter \n")
  }
  if(is.na(centerR)){ 
    stop("scaleOutlier: center parameter is NA \n")
  }
  if(is.na(scaleR)){ 
    stop("scaleOutlier: scale parameter is NA \n")
  }
  if(scaleR == 0){ 
    if(noScaleIf0){
      scaleR <- FALSE
    }else{
      stop("scaleOutlier: scale parameter is 0 \n")
    }
  }
  
  return( scale(x, center = centerR, scale = scaleR) )
  
}


#----------------------------------------------------------------------
### scaleOutlier.R ends here
