### qqplot2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: aug 30 2017 (09:26) 
## Version: 
## last-updated: maj  9 2018 (15:40) 
##           By: Brice Ozenne
##     Update #: 40
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * Documentation - qqplot2
#' @title QQplot for the Residuals of a LVM Object
#' @description QQplot for the residuals of a lvm object.
#' 
#' @name qqplot2
#'
#' @param object a lvm model.
#' @param residuals [matrix] the residuals relative to each endogenous variable.
#' @param variables the variable for which the residuals should be displayed.
#' @param name.model [character vector] character string to be displayed before the variable name in the title of the plot.#' 
#' @param mfrow how to divide the window. See \code{par}.
#' @param mar [numeric vector] the number of lines of margin
#' to be specified on the four sides of the plot (bottom, left, top, right).
#' @param type the function used to display the qqplot. Can be qqtest or qqnorm.
#' @param centralPercents argument passed to \code{qqtest}. See the help of \code{\link{qqtest}}.
#' @param ... additional arguments to be passed to qqtest.
#' 
#' @details 
#' Simulation is based on a multivariate truncated normal law (even though it is not satifying for the variance components)
#' 
#' @return a data frame/cvlvm object containing the convergence status (by default 0 indicates successful convergence, see ?optim), the value of the log-likelihood and the estimated parameters (in columns) for each initialization (in rows)
#' 
#' @examples
#' library(lava)
#' m <- lvm(list(y~v1+v2+v3+v4,c(v1,v2,v3,v4)~x))
#' latent(m) <- ~ x
#'
#' set.seed(10)
#' dd <- sim(m,100) ## Simulate 100 observations from model
#' e <- estimate(m, dd) ## Estimate parameters
#'
#' coef(e)
#' 
#' qqplot2(e)
#' 
#' @export
qqplot2 <- function (object, ...) {
  UseMethod("qqplot2", object)
}

## * qqplot2.lvmfit
#' @rdname qqplot2
#' @export
qqplot2.lvmfit <- function(object, variables = NULL, residuals = NULL, mfrow = NULL, mar = c(2,2,2,2),
                           type = "qqtest", name.model = "", centralPercents = 0.95,...){

    if(is.null(residuals)){
        residuals <- stats::predict(object, residual = TRUE)
    }
    name.vars <- colnames(residuals)
    
    if(!is.null(variables)){
        if(any(variables %in% name.vars == FALSE)){
            stop("unknown variable(s): ",paste(variables[variables %in% name.vars == FALSE], collapse = " "),"\n",
                 "endogenous variables: ",paste(endogenous(object), collapse = " "),"\n",
                 "latent variables: ",paste(latent(object), collapse = " "),"\n")
        }
        residuals <- residuals[,variables,drop=FALSE]
        name.vars <- variables
    }
    if(type %in% c("qqtest","qqnorm") == FALSE){
        stop("wrong specification of type \n",
             "must be \"qqtest\" or \"qqnorm\" \n")
    }

    n.var <- NCOL(residuals)       
    if(is.null(mfrow)){
        mfrow <- c(round(sqrt(n.var)), ceiling(n.var/round(sqrt(n.var))))
    }
    op <- graphics::par(mfrow = mfrow, mar = mar)
    sapply(1:n.var, function(row){
        resid <- stats::na.omit(residuals[,row])
        iMain <- paste0(name.model,": ",name.vars[row])
        if(all(resid < 1e-5)){
            graphics::plot(0,0, col = "white", axes = FALSE, xlab = "", ylab = "", main = iMain)
            graphics::text(0,0,"all residuals < 1e-5")
        }else if(type == "qqtest"){
            qqtest::qqtest(resid, main = iMain,
                           centralPercents = centralPercents,
                           ...)
        }else if(type == "qqnorm"){
            stats::qqnorm(residuals[,row], main = iMain)
        }
    })
    graphics::par(op)

    return(invisible(residuals))
}

## * qqplot2.multigroupfit
qqplot2.multigroupfit <- function(object, residuals = NULL, name.model = NULL, ...){
    if(is.null(residuals)){
        residuals <- stats::residuals(object)    
    }
    n.group <- length(residuals)
    if(is.null(name.model)){
        if(!is.null(object$model$names)){
            name.model <- object$model$names
        }else{
            name.model <- 1:n.group
        }
    }

    
    out <- vector(mode = "list", length = n.group)
    for(iG in 1:n.group){
        dev.new()
        out[[iG]] <- qqplot2.lvmfit(object,
                                    name.model = name.model[iG],
                                    residuals = residuals[[iG]],
                                    variables = colnames(residuals[[iG]]), ...)  
    }

    return(invisible(out))
}

##----------------------------------------------------------------------
### qqplot2.R ends here
