### plot.detectOutlier.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt  5 2017 (10:16) 
## Version: 
## last-updated: okt  5 2017 (10:16) 
##           By: Brice Ozenne
##     Update #: 3
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:


## * plot.detectOutlier
#' @title plot function for detectOutlier objects
#' @description plot function for detectOutlier objects.
#' 
#' @param x an object of class detectOutlier.
#' @param ... not used.
#'
#' @method plot detectOutlier
#' @export
plot.detectOutlier <- function(x, ...){
  
  if(length(x$display)>0){
    do.call(x$display$method, args = x$display$args)
  }else{
    message("no available plot \n")
  }
  
}


#----------------------------------------------------------------------
### plot.detectOutlier.R ends here
