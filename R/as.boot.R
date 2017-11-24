### as.boot.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 22 2017 (09:41) 
## Version: 
## Last-Updated: nov 22 2017 (17:49) 
##           By: Brice Ozenne
##     Update #: 15
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
  function(object) UseMethod("as.boot")

## * as.boot.bootReg
#' @rdname as.boot
#' @export
as.boot.bootReg <- function(object){

    n.boot.effective <- NROW(object$boot.estimate)
    n.data <- NROW(object$data)
    out <- list(t0 = object$estimate,
                t = object$boot.estimate,
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