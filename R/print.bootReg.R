### print.bootReg.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 22 2017 (15:44) 
## Version: 
## Last-Updated: nov 22 2017 (16:18) 
##           By: Brice Ozenne
##     Update #: 6
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * print.bootReg
#' @method print bootReg
#' @export
print.bootReg <- function(x, ...){

    out <- summary(x, p.value = FALSE)
    return(invisible(x))
    
}
##----------------------------------------------------------------------
### print.bootReg.R ends here
