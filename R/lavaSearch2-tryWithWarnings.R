### tryWithWarnings.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: aug 29 2017 (09:52) 
## Version: 
## last-updated: nov 21 2017 (19:29) 
##           By: Brice Ozenne
##     Update #: 18
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:


##' @title Run an expression and catch warnings and errors
##' @description Similar to \code{try} but also returns warnings.
##' 
##' @param expr the line of code to be evaluated
##' @details
##' from https://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function
##' 
##' @return a list containing:
##' \itemize{
##' \item value the result of the evaluation of the expression
##' \item warnings warning(s) generated during the evaluation of the expression
##' \item error error generated during the evaluation of the expression
##' }
##'
##' @examples
##' FctTest <- function(x){
##'   return(log(x))
##' }
##' tryWithWarnings(FctTest(-1))
##' tryWithWarnings(FctTest(1))
##' tryWithWarnings(FctTest(xxxx))
##'' '
##' @author Brice Ozenne
##' @export
tryWithWarnings <- function(expr) {
    myWarnings <- NULL
    myError <- NULL
    FCT.envir <- environment()
    wHandler <- function(w) {
        assign("myWarnings", value = c(myWarnings, list(w)), envir = FCT.envir)
        invokeRestart("muffleWarning")
    }
    eHandler <- function(e) {
        assign("myError", value = e, envir = FCT.envir)
        return(NULL)
    }
    val <- withCallingHandlers(tryCatch(expr, error = eHandler), warning = wHandler)
    
    out <- list(value = val,
                warnings = myWarnings,
                error = myError)
    
  return(out)
}


#----------------------------------------------------------------------
### tryWithWarnings.R ends here
