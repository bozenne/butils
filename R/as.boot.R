### as.boot.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 22 2017 (09:41) 
## Version: 
## Last-Updated: apr 11 2018 (15:31) 
##           By: Brice Ozenne
##     Update #: 26
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:
#' @title Convertion to boot object 
#' @description Convertion to boot object
#' @name as.boot
#'  
#' @param object a bootReg object
#' @param index [optional] index of the observation to keep.
#' 
#' @examples
#' #### data  ####
#' n <- 1e2
#' set.seed(10)
#' df.data <- data.frame(Y = rnorm(n),
#'                      group = gl(3, 5, n, labels = c("Ctl","Trt","Neu")),
#'                      gender = gl(2, 5, n, labels = c("Female","Male"))[sample.int(n)]
#'                      )
#' 
#' #### lm ####
#' m.lm <- lm(Y ~ group*gender, data = df.data)
#' resBoot <- bootReg(m.lm, n.boot = 1e1)
#' as.boot(resBoot)
#' 
#' @export
`as.boot` <-
  function(object, index) UseMethod("as.boot")

## * as.boot.bootReg
#' @rdname as.boot
#' @export
as.boot.bootReg <- function(object, index = NULL){

    if(is.null(index)){
        index <- 1:length(object$estimate)
    }
    
    n.boot.effective <- NROW(object$boot.estimate)
    n.data <- NROW(object$data)
    out <- list(t0 = object$estimate[index],
                var.t0 = object$stdError[index]^2,
                t = object$boot.estimate[,index,drop=FALSE],
                var.t = object$boot.stdError[,index,drop=FALSE]^2,
                R = n.boot.effective,
                data = object$data,
                seed = object$.Random.seed,
                statistic = object$FUN.estimate,
                sim = "ordinary",
                call = quote(boot(data = XX, statistic = XX, R = XX)),
                stype = "i",
                strata = object$strata,
                weights = rep(1/n.data, n.data)
                )
    class(out) <- "boot"
    return(out)
 
}
##----------------------------------------------------------------------
### as.boot.R ends here
