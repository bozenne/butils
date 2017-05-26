##' Compute partial residuals for linear and linear mixed models.
##'
##' @title Plot regression lines
##' @param model Model object (e.g. \code{lm})
##' @param var Variable relative to which the partial residual will be displayed.
##' @param conditional are the predictions conditional to the random effect? (if any)
##' @param interval Type of interval calculation ("confidence" or "prediction").
##' @param level Level of confidence limits (default 95\%)
##' @param npoints Length of the vector of unique values relative to each continuous variable.
##' @param quantile.range the quantiles (for the continous covariates) between which the fitted values will be computed.
##' @param FUN.predict, Optional predict-function used to calculate confidence limits and predictions (see details).
##' @param FUN.df Optional function returning the residual degree of freedoms.
##' @param \dots additional arguments to lower level functions
##' 
##' @details
##' When using mixed models, the confidence and prediction intervals ignore the uncertainty of the covariance parameters.
##' To avoid that the user can set the argument \code{predictfun} to the appropriate function to compute prediction and confidence intervals.
##' This function should take the following arguments: object, newdata, level, interval, se.fit (will be automatically be set to TRUE).
##' It should output a list with an element called fit that contains a matrix with one row for each line of newdata
##' and three columns:
##' - \code{fit} containing the punctual estimate
##' - \code{upr} containing the upper limit of the confidence interval
##' - \code{lwr} containing the lower limit of the confidence interval
##' 
##' @return list with following members:
##' \item{data}{Original dataset with an additional column containing the partial residuals}
##' \item{partialFit}{Fitted values according to the range of values of var}
##' \item{var}{variable(s) for which the partial residuals have been computed}
##' \item{level}{Level of confidence}
##' \item{name.Y}{name of the response variable}
##' @author Brice Ozenne
##' @seealso \code{plot.partialResiduals} for a graphical display and more examples.
##' @export
##' @examples
##' library(lava)
##' set.seed(10)
##' m.lvm <- lvm(Y~X1+X2+Id+Id2)
##' categorical(m.lvm, K = 5) <- ~Id+Id2
##' d <- sim(n = 1e2, m.lvm)
##'  
##' m <- lm(Y~X1+X2, data = d)
##' pres1 <- calcPartialResiduals(m, var = "X1")
##' pres2 <- calcPartialResiduals(m, var = c("X1","X2"))
##'
##' 
##' @keywords hplot, regression
calcPartialResiduals <- function(model,var,
                                 conditional = FALSE,
                                 interval = "confidence",                                 
                                 level = 0.95,
                                 npoints = 100,
                                 quantile.range = c(0,1),                                 
                                 FUN.df,
                                 FUN.predict,
                                 ...) {

    model.class <- class(model)
    attributes(model.class) <- NULL
    if (!identical(model.class,"lmerMod") && !identical(model.class,"lme") && !identical(model.class,"gls") && !identical(model.class,"lm")){
        stop("calcPartialResiduals can only deal with lm, lme, and lmerMod objects \n")
    }

    out <- list()

    # {{{ define class specific functions
    if(model.class %in% c("lme","lmerMod")){
        FUN.coef <- function(object, ...){
            fixef(object, ...)
        }
    }else if(model.class %in% c("lm","gls")){
        FUN.coef <- function(object, ...){
            coef(object, ...)
        }
    }

    if(model.class == "lme"){
        FUN.ranef <- function(object, ...){
            ranef(object, ...)
        }
    }else if(model.class == "lmerMod"){
        FUN.ranef <- function(object, ...){
            ranef(object, ...)
        }
    }else if(model.class %in% c("lm","gls","lme")){
        FUN.ranef <- function(object, ...){
            coef(object, ...)
        }
    }
    
    
    if(model.class=="lmerMod"){
        FUN.vcov <- function(object, ...){
            as.matrix(stats::vcov(object, ...))  
        }
    }else if(model.class %in% c("lm","gls","lme")){
        FUN.vcov <- function(object, ...){
            stats::vcov(object, ...)
        }
    }

    FUN.sigma2 <- function(object, ...){
        summary(object, ...)$sigma^2
    }

    FUN.model.matrix <- function(model, data){
        model.matrix(model, data = data, fixed.only = FALSE)
    }

    if(model.class %in% c("lme","gls")){
        FUN.model.frame <- function(model, ...){
            nlme::getData(model)
        }
    }else{
        FUN.model.frame <- function(model, ...){
            model.frame(model, fixed.only = FALSE)
        }
    }

    if(missing(FUN.df)){
        if(model.class %in% c("lmerMod","lme","gls")){
            FUN.df <- function(model, level){
                NULL
            }
        }else if(model.class=="lm"){
            FUN.df <- function(model, level){
                model$df.residual
            }

        }
    }
    
    FUN.formula <- function(model){
        formula(terms(model))
    }
    # }}}

    # {{{ normalize input

    interval <- match.arg(interval, c("confidence","prediction"))

    design.df <- FUN.model.frame(model)
    
    design.numeric <- sapply(design.df, is.numeric)
    design.factor <- sapply(design.df, is.factor)

    design.NonNumeric <- names(which(design.numeric==FALSE))
    if(length(design.NonNumeric)>0){
        xlevels <- sapply(design.NonNumeric, function(x){
            list(levels(as.factor(design.df[[x]])))
        })
    }
    
    name.Y <- all.vars(formula(model))[1]
    name.X.df <- setdiff(colnames(design.df),name.Y)
    if (any(var %in% name.X.df == FALSE) ) {
        stop("Unknown variable(s): \"",paste(var[var %in% name.X.df == FALSE], collapse = "\" \""),"\" \n",
             "argument \'var\' should be in \"",paste(name.X.df,collapse = "\" \""),"\" \n")
    }

    # }}}
  
    # {{{ additional tests
    if (inherits(model,"lmerMod") && missing(FUN.predict)) {
        warning(interval," intervals may not be reliable (see ?lme4:::predict.merMod) \n",
                "Consider specifying predictfun (e.g. using lme4::bootMer) \n")
    }

    
    # }}}

    # {{{ Reference level
    n <- NROW(design.df)
    p <- length(name.X.df)
    newdata.fit <- as.data.frame(design.df)

    var.set0 <- names(design.df)[which(names(newdata.fit) %in% var == TRUE)]    
    if(length(var.set0)){
        for(i0 in var.set0){ # i0 <- var.set0[1]
            if(design.numeric[i0]){
                newdata.fit[[i0]][] <- 0
            }else{
                newdata.fit[[i0]][] <- xlevels[[i0]][1]
            }
            
        }
    }
    # }}}

    # {{{ intercept
    # if intercept always remove it from the fitted value
    # if no intercept:
    # - no factor variable: 0 interceot
    # - first factor variable: not of interest treat it as a standard intercept
    #                              of interest remove it only for partial residuals
    model.formula <- FUN.formula(model)
    beta <- FUN.coef(model)
    if(attr(terms(model.formula),"intercept")==1){
        intercept.data <- 0       
        intercept.fit <- beta["(Intercept)"]        
    }else{
        if(any(design.factor)){ # find reference level
            design.mat <- FUN.model.matrix(model.formula, data = design.df[1,,drop=FALSE])
            attr.assign <- attr(design.mat, "assign")
            indexRef <- which(tapply(attr.assign, attr.assign, function(x){any(duplicated(x))}))
            if(attr(terms(model.formula), "term.labels")[indexRef] %in% var){
                intercept.data <- beta[match(indexRef, attr.assign)]
                intercept.fit <- 0
            }else{
                intercept.data <- 0
                intercept.fit <- beta[match(indexRef, attr.assign)]
            }
        }else{
            intercept.data <- 0
            intercept.fit <- 0
        }
    }
    # }}}

    # {{{ partial residuals
    if(!missing(FUN.predict)){
        design.df$pFit <- do.call(FUN.predict, args = list(model, newdata = newdata.fit, conditional = conditional, type = "response"))$fit[,"fit"]
    }else{        
        design.mat <- model.matrix(model.formula, data = newdata.fit)
        design.df$pFit <- as.numeric(design.mat %*% beta)
        if(conditional){ # https://stackoverflow.com/questions/25538199/design-matrix-for-mlm-from-librarylme4-with-fixed-and-random-effects
            stop("no implemented - specify FUN.predict \n")
        }
    }
    design.df$pResiduals <- newdata.fit[[name.Y]] - (design.df$pFit-intercept.data)
    # }}}
    
    # {{{ partial fitted values
    ## newdata
    ls.forGrid <- lapply(name.X.df, function(x){
        if(x %in% var == FALSE){
            if(design.numeric[x]){
                return(0)
            }else {
                return(xlevels[[x]][1])
            }
        }else if(design.numeric[x]){
            return(seq(quantile(design.df[[x]], na.rm=TRUE, probs = quantile.range[1]),
                       quantile(design.df[[x]], na.rm=TRUE, probs = quantile.range[2]),
                       length.out = npoints))
        }else{ 
            return(unique(design.df[[x]]))
        }
    })

    grid.predict <- expand.grid(ls.forGrid, KEEP.OUT.ATTRS = FALSE)
    names(grid.predict) <- name.X.df
    if(any(design.factor)){         # restaure factor
        for(iVar in names(which(design.factor))){
            grid.predict[[iVar]] <- factor(grid.predict[[iVar]], levels = xlevels[[iVar]])
        }
    }
    
    if(!missing(FUN.predict)){
        res.predict <- do.call(FUN.predict, args = list(model,
                                                        se.fit = TRUE,
                                                        conditional = conditional,
                                                        newdata = grid.predict,
                                                        interval = interval,
                                                        level = level,
                                                        ...))

        grid.predict$fit <- res.predict$fit[,"fit"]
        grid.predict$fit.lower <- res.predict$fit[,"lwr"]
        grid.predict$fit.upper <- res.predict$fit[,"upr"]
    }else{
        if(conditional){ # https://stackoverflow.com/questions/25538199/design-matrix-for-mlm-from-librarylme4-with-fixed-and-random-effects
            stop("no implemented - specify FUN.predict \n")
        }
        
        grid.predict[[name.Y]] <- 0
        design.grid.predict <- model.matrix(model.formula, grid.predict)
        grid.predict$fit <- as.numeric(design.grid.predict %*% beta)

        model.df <- FUN.df(model,level)       
        if(is.null(model.df)){
            z <- qnorm(1 - (1 - level)/2)
        }else{
            z <- qt(1 - (1 - level)/2, df = model.df)
        }
        
        model.vcov <- FUN.vcov(model)
        if(interval == "prediction"){
            model.sigma2 <- FUN.sigma2(model)
        }else{
            model.sigma2 <- 0
        }
        
        se.tempo <- apply(design.grid.predict, 1, function(x){
            sqrt(rbind(x) %*% model.vcov %*% cbind(x) + model.sigma2)
        })

        grid.predict$fit.lower <- grid.predict$fit - z * se.tempo
        grid.predict$fit.upper <- grid.predict$fit + z * se.tempo        
    }
    grid.predict$fit <- grid.predict$fit - intercept.fit
    grid.predict$fit.lower <- grid.predict$fit.lower - intercept.fit
    grid.predict$fit.upper <- grid.predict$fit.upper - intercept.fit

    # }}}

    out <- list(data = as.data.table(design.df),                
                partialFit = as.data.table(grid.predict),
                var = var,
                interval = interval,
                level = level,
                name.Y = name.Y)
    class(out) <- "partialResiduals"
    return(out)
}


# {{{ prediction functions
#' @title prediction function for calcPartialResiduals
#' @rdname FUN.predict
#' @description Suggestion of functions to obtain confidence intervals with mixed models
#'
#' @param object a model
#' @param newdata the dataset used to make the predictions
#' @param level the confidence level for the intervals
#' @param interval the type of interval, prediction or confidence
#' @param conditional are the predictions conditional to the random effect? (if any)
#' @param ... not used
#' 
#' @details
#' predict_merTools is based on the \code{predictInterval} function of the merTools package
#' and is only compatible with \code{lmer} models.
#'
#' predict_AICcmodavg is based on the \code{predictSE} function of the AICcmodavg package
#' and is compatible with \code{lmer}, \code{lme}, and \code{gls} models. It can only be used to obtain confidence intervals.
#'
#' See the details of the documentation of these functions for more information about the validity of the intervals.
#' 

# {{{ predict_merTools
#' @rdname FUN.predict
#' @export
predict_merTools <- function(object, newdata, level,conditional, interval, ...){
    if(missing(interval)){interval <- FALSE}
    if(missing(level)){level <- 0.95}
    res <- merTools::predictInterval(object, newdata = newdata, level = level,
                                     type = "linear.prediction", n.sims = 1e3,
                                     which = if(conditional){"full"}else{"fixed"},
                                     include.resid.var = (interval=="prediction"))
    return(list(fit = res))
}
# }}}

# {{{ predict_AICcmodavg 
#' @rdname FUN.predict
#' @export
predict_AICcmodavg <- function(object, newdata, level, conditional, interval, ...){
    if(!missing(interval) && interval == "prediction"){
        stop("no implemented for prediction intervals \n")
    }    
    res <- AICcmodavg::predictSE(object, newdata = newdata,
                                 type = "response", level = as.numeric(conditional))
    if(missing(level)){
        M.res <- cbind(fit = res$fit)
    }else{
        z <- qnorm(1 - (1 - level)/2)
        M.res <- cbind(fit = res$fit,
                       lwr = res$fit - z * res$se.fit,
                       upr = res$fit + z * res$se.fit)

    }
    return(list(fit = M.res))
}
# }}}

# }}}
