## * documentation - calcPartialResiduals
##' @title Compute partial residuals for linear and linear mixed models.
##' @description Compute partial residuals for linear and linear mixed models.
##'
##' @name calcPartialResiduals
##' 
##' @param model Model object (e.g. \code{lm})
##' @param var Variable relative to which the partial residual will be displayed.
##' @param keep.intercept should the partial residual be computed keeping the contribution of the reference level?
##' @param conditional are the predictions conditional to the random effect? (if any)
##' @param interval Type of interval calculation ("confidence" or "prediction").
##' @param level Level of confidence limits (default 95\%)
##' @param npoints Length of the vector of unique values relative to each continuous variable.
##' @param quantile.range the quantiles (for the continous covariates) between which the fitted values will be computed.
##' @param FUN.df Optional function returning the residual degree of freedoms.
##' @param \dots additional arguments to lower level functions
##' 
##' @details
##' In a linear model:
##' \deqn{Y = \alpha + \beta X + \gamma Z + \varepsilon}
##' The partial residuals are defined by:
##' \deqn{\varepsilon_{X} = \beta X + \varepsilon}
##' or
##' \deqn{\varepsilon_{X} = \alpha + \beta X + \varepsilon}
##' depending on the value of the argument \code{keep.intercept}.
##' 
##' When using mixed models, the confidence and prediction intervals
##' ignore the uncertainty of the covariance parameters.
##' 
##' @return list with following members:
##' \item{data}{Original dataset with an additional column containing the partial residuals}
##' \item{partialFit}{Fitted values according to the range of values of var}
##' \item{var}{variable(s) for which the partial residuals have been computed}
##' \item{level}{Level of confidence}
##' \item{name.Y}{name of the response variable}
##'
##' @seealso \code{plot.partialResiduals} for a graphical display and more examples.'
##' 
##' @author Brice Ozenne
##' @examples
##' library(lava)
##' set.seed(10)
##' m.lvm <- lvm(Y~X1+X2+Id+Id2)
##' categorical(m.lvm, K = 5) <- ~Id+Id2
##' d <- lava::sim(n = 1e2, m.lvm)
##'  
##' m <- lm(Y~X1+X2, data = d)
##' pres1 <- calcPartialResiduals(m, var = "X1")
##' pres2 <- calcPartialResiduals(m, var = c("X1","X2"))
##'
##' 
##' @keywords regression

## * function - calcPartialResiduals
##' @rdname calcPartialResiduals
##' @export
calcPartialResiduals <- function(model,var,
                                 keep.intercept = FALSE,
                                 conditional = FALSE,
                                 interval = "confidence",                                 
                                 level = 0.95,
                                 npoints = 100,
                                 quantile.range = c(0,1),                                 
                                 FUN.df,
                                 ...) {

    model.class <- class(model)
    attributes(model.class) <- NULL
    if (!identical(model.class,"lmerMod") && !identical(model.class,"lme") && !identical(model.class,"gls") && !identical(model.class,"lm")){
        stop("calcPartialResiduals can only deal with lm, lme, and lmerMod objects \n")
    }

    out <- list()
#
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
      ff <- try(formula(terms(model)), silent = TRUE)
      if("try-error" %in% class(ff)){
        ff <- evalInParentEnv(model$call[[2]])
      }
      
      if("formula" %in% class(ff) == FALSE){
       stop("Unable to extract the formula from the model \n")
      }
      return(ff)
    }
    # }}}

    # {{{ normalize input

    interval <- match.arg(interval, c("confidence","prediction"))

    design.df <- as.data.table(extractData(model, design.matrix = FALSE))
    
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
    ## if (inherits(model,"lmerMod")) {
    ##     warning(interval," intervals may not be reliable (see ?lme4:::predict.merMod) \n")
    ## }

    if(conditional){ # https://stackoverflow.com/questions/25538199/design-matrix-for-mlm-from-librarylme4-with-fixed-and-random-effects
        stop("no yet implemented \n")
    }

    # }}}

    # {{{ Reference level
    n <- NROW(design.df)
    p <- length(name.X.df)
    newdata.fit <- as.data.frame(design.df)

    var.set0 <- names(newdata.fit)[which(names(newdata.fit) %in% var == TRUE)]
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
    model.formula <- FUN.formula(model)
    beta <- FUN.coef(model)
    test.intercept <- as.logical(attr(terms(model.formula),"intercept"))
    
    if(test.intercept){
      name.intercept <- "(Intercept)"
    }else if(any(design.factor)){ # find reference level
            design.mat <- FUN.model.matrix(model.formula, data = design.df[1,,drop=FALSE])
            attr.assign <- attr(design.mat, "assign")
            indexRef <- which(tapply(attr.assign, attr.assign, function(x){any(duplicated(x))}))
            name.intercept <- names(beta)[match(indexRef, attr.assign)]
    }else{
      name.intercept <- NULL
    }
    # }}}

    # {{{ partial residuals
    design.mat <- model.matrix(model.formula, data = newdata.fit)    
    if(!is.null(name.intercept) && keep.intercept){
      # add the residual due to the intercept
      # since it is among the variable of interest
      design.mat[,name.intercept] <- 0
    }
    design.df$pFit <- as.numeric(design.mat %*% beta)
    design.df$pResiduals <- newdata.fit[[name.Y]] - design.df$pFit
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
            return(seq(stats::quantile(design.df[[x]], na.rm=TRUE, probs = quantile.range[1]),
                       stats::quantile(design.df[[x]], na.rm=TRUE, probs = quantile.range[2]),
                       length.out = npoints))
        }else{ 
            return(unique(design.df[[x]]))
        }
    })

    ## create the data.frame
    grid.predict <- expand.grid(ls.forGrid, KEEP.OUT.ATTRS = FALSE)
    names(grid.predict) <- name.X.df
    grid.predict[[name.Y]] <- 0
    ## restaure factor variables
    if(any(design.factor)){         
        for(iVar in names(which(design.factor))){
            grid.predict[[iVar]] <- factor(grid.predict[[iVar]], levels = xlevels[[iVar]])
        }
    }
    ## convert to design matrix
    design.grid.predict <- model.matrix(model.formula, grid.predict)
    ## remove intercept
    if(!is.null(name.intercept) && keep.intercept == FALSE){
        # remive the intercept from the prediction
        # when it is not among the variable of interest
        design.grid.predict[,name.intercept] <- 0
    }
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


## * predict functions
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

