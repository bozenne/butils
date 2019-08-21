### qqplot2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: aug 30 2017 (09:26) 
## Version: 
## last-updated: jun 27 2019 (09:58) 
##           By: Brice Ozenne
##     Update #: 108
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
#' @param name.model [character vector] character string to be displayed before the variable name in the title of the plot.
#' If \code{NULL} and \code{object} is a list, the names of the list will be used to define the title of the plots.
#' @param mfrow how to divide the window. See \code{par}.
#' @param mar [numeric vector] the number of lines of margin
#' to be specified on the four sides of the plot (bottom, left, top, right).
#' @param plot [logical] should the graphic be displayed?
#' @param qq.type the function used to display the qqplot. Can be qqtest or qqnorm.
#' @param centralPercents argument passed to \code{qqtest}. See the help of \code{qqtest}.
#' @param ... additional arguments to be passed to qqtest.
#' 
#' @details 
#' Simulation is based on a multivariate truncated normal law (even though it is not satifying for the variance components)
#' 
#' @return a data frame/cvlvm object containing the convergence status (by default 0 indicates successful convergence, see ?optim), the value of the log-likelihood and the estimated parameters (in columns) for each initialization (in rows)
#' 
#' @examples
#' #### lvm object ####
#' library(lava)
#' m <- lvm(list(y~v1+v2+v3+v4,c(v1,v2,v3,v4)~x))
#' latent(m) <- ~ x
#'
#' set.seed(10)
#' dd <- sim(m,100) ## Simulate 100 observations from model
#' e <- estimate(m, dd) ## Estimate parameters
#'
#' qqplot2(e)
#'
#' #### multigroup object ####
#' e2 <- estimate(list(m1=m,m2=m), split(dd, dd$v1>0)) ## Estimate parameters
#'
#' qqplot2(e2)
#' 
#' @export
qqplot2 <- function (object, ...) {
  UseMethod("qqplot2", object)
}

## * qqplot2.lvmfit
#' @rdname qqplot2
#' @export
qqplot2.lvmfit <- function(object, variables = NULL, residuals = NULL,
                           plot = TRUE, mfrow = NULL, mar = c(2,2,2,2),
                           qq.type = "qqtest", name.model = "", centralPercents = 0.95, ...){
    
    ## ** prepare
    if(is.null(residuals)){
        residuals <- stats::predict(object, residual = TRUE)
    }
    if(is.null(variables)){
        variables <- colnames(residuals)
    }
    
    out <- .qqplot2(variables = variables,
                    residuals = residuals,
                    plot = plot,
                    mfrow = mfrow,
                    mar = mar,
                    qq.type = qq.type,
                    name.model = name.model,
                    centralPercents = centralPercents)
    
    return(invisible(out))
}

## * qqplot2.list
#' @rdname qqplot2
#' @export
qqplot2.list <- function(object, 
                         plot = TRUE, mfrow = NULL, mar = c(2,2,2,2),
                         qq.type = "qqtest", name.model = "", centralPercents = 0.95, ...){

    n.list <- length(object)
    ls.residuals <- vector(mode = "list", length = n.list)
    variables <- vector(mode = "character", length = n.list)
    vec.n <- vector(mode = "numeric", length = n.list)

    for(iL in 1:n.list){
        ls.residuals[[iL]] <- stats::residuals(object[[iL]], ...)
        variables[iL] <- all.vars(formula(object[[iL]]))[1]
        vec.n[iL] <- length(ls.residuals[[iL]])
    }
    if(any(vec.n != vec.n[iL])){
        stop("The number of residuals differ between models \n")
    }
    if(any(duplicated(variables))){
        if(length(unique(name.model))!=length(variables) || any(duplicated(name.model))){
            stop("Non unique response variable names \n")
        }else{
            variables <- name.model
            name.model <- ""
        }
    }
    
    M.residuals <- do.call("cbind",ls.residuals)
    colnames(M.residuals) <- variables

    out <- .qqplot2(variables = variables,
                    residuals = M.residuals,
                    plot = plot,
                    mfrow = mfrow,
                    mar = mar,
                    qq.type = qq.type,
                    name.model = name.model,
                    centralPercents = centralPercents)
    return(invisible(out))
    
}

## * qqplot2.multigroupfit
#' @rdname qqplot2
#' @export
qqplot2.multigroupfit <- function(object, residuals = NULL, name.model = NULL, plot = TRUE, ...){
    if(is.null(residuals)){
        residuals <- stats::residuals(object)    
    }
    n.group <- length(residuals)
    if(is.null(name.model)){
        if(!is.null(object$model$names)){
            name.model <- paste0(object$model$names,": ")
        }else{
            name.model <- paste0("model ",1:n.group,": ")
        }
    }

    
    out <- vector(mode = "list", length = n.group)
    for(iG in 1:n.group){
        if(plot == TRUE){ grDevices::dev.new() }
        out[[iG]] <- qqplot2.lvmfit(object,
                                    name.model = name.model[iG],
                                    residuals = residuals[[iG]],
                                    variables = colnames(residuals[[iG]]),
                                    plot = plot,
                                    ...)  
    }

    return(invisible(out))
}

## * .qqplot2
.qqplot2 <- function(variables, residuals,
                     plot, mfrow, mar,
                     qq.type, name.model, centralPercents, ...){

    ## ** check arguments
    if(any(variables %in% colnames(residuals) == FALSE)){
        txt <- paste(variables[variables %in% colnames(residuals) == FALSE], collapse = "\" \"")
        stop("unknown variable(s): \"",,"\"\n")
    }
    if(qq.type %in% c("qqtest","qqnorm") == FALSE){
        stop("wrong specification of qq.type \n",
             "must be \"qqtest\" or \"qqnorm\" \n")
    }

    ## ** initialize
    residuals <- as.data.frame(residuals)[,variables,drop=FALSE]
    n.var <- NCOL(residuals)       
    if(is.null(mfrow)){
        mfrow <- c(round(sqrt(n.var)), ceiling(n.var/round(sqrt(n.var))))
    }

    ## ** wraper
    warper_display <- function(envir){
        op <- graphics::par(mfrow = envir$mfrow, mar = envir$mar)
        sapply(1:envir$n.var, function(row){
            resid <- stats::na.omit(envir$residuals[,row])
            if(length(envir$name.model)>1){
                iMain <- envir$name.model[row]
            }else{
                iMain <- paste0(envir$name.model,envir$variables[row])
            }
            if(all(resid < 1e-5)){
                graphics::plot(0,0, col = "white", axes = FALSE, xlab = "", ylab = "", main = iMain)
                graphics::text(0,0,"all residuals < 1e-5")
            }else if(envir$qq.type == "qqtest"){
                requireNamespace("qqtest")
                qqtest::qqtest(resid, main = iMain,
                               centralPercents = envir$centralPercents,
                               ...)
            }else if(envir$qq.type == "qqnorm"){
                stats::qqnorm(resid, main = iMain)
            }
        })
        graphics::par(op)
    }

    ## ** gather
    ls.data <- list(mfrow = mfrow,
                    mar = mar,
                    n.var = n.var,
                    residuals = residuals,
                    name.model = name.model,
                    variables = variables,
                    centralPercents = centralPercents,
                    qq.type = qq.type)
    
    if(plot){
        warper_display(ls.data)
    }
    
    ## export
    return(invisible(list(plot = warper_display,
                          data = ls.data)))
}


##----------------------------------------------------------------------
### qqplot2.R ends here
