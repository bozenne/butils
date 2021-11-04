## * documentation - partialResiduals
##' @title Compute partial residuals for linear and linear mixed models.
##' @description Compute partial residuals for linear and linear mixed models.
##' @name partialResiduals
##' 
##' @param model Model object (e.g. \code{lm})
##' @param var [character] Independent variable(s) whose effect should be kept in the partial residuals.
##' @param keep.intercept [logical] should the partial residual be computed keeping the contribution of the reference level?
##' Can also a vector of lenght 2, the second element indicating whether the uncertainty calculation should account for the estimation of the intercept.
##' @param conditional [logical] should the partial residuals be computed without the random effects? (if any)
##' @param interval [character] Type of interval calculation ("confidence" or "prediction").
##' @param level [numeric, 0-1] Level of confidence limits (default 95\%)
##' @param npoints [integer] Length of the vector of unique values relative to each continuous variable.
##' @param quantile.range [numeric vector, 0-1] the quantiles (for the continous covariates) between which the fitted values will be computed.
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
##' The X matrix is defined contains the variables defined by the \code{var} argument.
##'
##' Confidence intervals are only valid in homoschedastic models.
##' When using mixed models, the confidence and prediction intervals
##' ignore the uncertainty of the covariance parameters/random effects.
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
##' m.lvm <- lvm(Y~X1+X2+Id)
##' categorical(m.lvm, K = 5) <- ~Id
##' d <- lava::sim(n = 1e2, m.lvm)
##' d$Id <- as.factor(d$Id)
##' m <- lm(Y~X1+X2+Id, data = d)
##'
##' #### 1- partial residuals relative to no variable matches the residuals
##' pres0 <- partialResiduals(m, var = NULL)
##' range(residuals(m) - pres0$pResiduals)
##'
##' #### 2- partial residuals regarding the variable X1
##' pres1 <- partialResiduals(m, var = "X1")
##' fit1 <- coef(m)["(Intercept)"] + coef(m)["X2"] * d$X2 +  coef(m)["Id"] * d$Id
##' range((d$Y - fit1) - pres1$pResiduals)
##' cor(d$X1, pres1$pResiduals)
##' lava::partialcor(~X2+Id,d[,c("Y","X1","X2","Id")])
##' 
##' \dontrun{
##' ## graphical display
##' autoplot(partialResiduals(m, var = "X1", keep.intercept = TRUE))
##' autoplot(partialResiduals(m, var = "X1", keep.intercept = FALSE))
##' autoplot(partialResiduals(m, var = "X1", keep.intercept = c(TRUE,FALSE)))
##' autoplot(partialResiduals(m, var = "Id", keep.intercept = TRUE))
##' }
##' 
##' #### 3-  partial residuals regarding both X1 and X2
##' pres2 <- partialResiduals(m, var = c("X1","X2"))
##' fit2 <- coef(m)["(Intercept)"] + d$Id * coef(m)["Id"]
##' range((d$Y - fit2) - pres2$pResiduals)
##'
##' \dontrun{
##' ## graphical display
##' autoplot(partialResiduals(m, var = c("X1","Id"), keep.intercept = TRUE))
##' autoplot(partialResiduals(m, var = c("X1","Id"), keep.intercept = FALSE))
##' autoplot(partialResiduals(m, var = c("X1","Id"), keep.intercept = c(TRUE,FALSE)))
##' }
##'
##' ## partial residuals binary variable
##' pres3 <- partialResiduals(m, var = "Id")
##' 
##' if(require(ggplot2)){
##'    pres3$Id.char <- factor(pres3$Id,
##'                                 levels = 0:4,
##'                                 labels = c("a","b","c","d","e"))
##'    gg <- ggplot(pres3, aes(y = pResiduals, group = Id.char, x = Id.char))
##'    gg <- gg + geom_boxplot()
##'    gg <- gg + geom_dotplot(binaxis = "y", stackdir = "center")
##'    gg
##' }
##'
##' ## partial residuals in presence of interactions
##' m.I <- lm(Y~X1*X2+Id, data = d)
##' pres.I <- partialResiduals(m.I, var = c("X1","X2"))
##' fit.I <- coef(m.I)["(Intercept)"] + d$Id * coef(m.I)["Id"]
##' range((d$Y - fit.I) - pres.I$pResiduals)
##'
##' 
##' 
##' @keywords regression

## * function - partialResiduals
##' @rdname partialResiduals
##' @export
partialResiduals <- function(model,var,
                             keep.intercept = FALSE,
                             conditional = FALSE,
                             interval = "confidence",                                 
                             level = 0.95,
                             npoints = 100,
                             quantile.range = c(0,1),                                 
                             FUN.df,
                             ...) {

    
    test.lm <- inherits(model,"lm")
    test.gls <- inherits(model,"gls")
    test.gam <- inherits(model,"gam")
    test.lme <- inherits(model,"lme")
    test.lmer <- inherits(model,"lmerMod")
    if(test.lm==FALSE && test.gls==FALSE && test.gam==FALSE && test.lme==FALSE && test.lmer==FALSE){
        stop("partialResiduals can only deal with lm, gls, gam, lme, and lmerMod objects \n")
    }
    
    out <- list()
    ## ** define class specific functions
    if(test.lme ||test.lmer){
        FUN.coef <- function(object, ...){
            fixef(object, ...)
        }
    }else if(test.lm || test.gam || test.gls){
        FUN.coef <- function(object, ...){
            coef(object, ...)
        }
    }

    if(conditional){
        if(test.lmer){
            FUN.ranef <- function(object, ...){
                predict(model, re.form = NULL, random.only = TRUE, type = "response")
            }
        }else if(test.lme){
            FUN.ranef <- function(object, ...){
                predict(model) - predict(model, level = 0)
            }
        }else if(test.lm || test.gam || test.gls){
            stop("Argument \'conditional=TRUE\' not applicable to lm or gls objects. \n")
        }
    }
    
    
    if(test.lmer){
        FUN.vcov <- function(object, ...){
            as.matrix(stats::vcov(object, ...))  
        }
    }else if(test.lm || test.gls || test.gam || test.lme){
        FUN.vcov <- function(object, ...){
            stats::vcov(object, ...)
        }
    }

    FUN.sigma2 <- function(object, ...){
        summary(object, ...)$sigma^2
    }

    if(test.gam){
        FUN.model.matrix <- function(model, formula, data){
            model.matrix(model, newdata = data, fixed.only = FALSE)
        }
    }else if(test.lmer){
        FUN.model.matrix <- function(model, formula, data){
            X <- model.matrix(formula, data = data, fixed.only = FALSE)

            ## remove extra column for random effect
            ff.random <- lme4::findbars(formula)
            query <- paste(paste0("^",unlist(lapply(ff.random, deparse))),collapse = "|")
            position.query <- grep(query, colnames(X))
            
            XX <- X[,-position.query,drop=FALSE]
            attr(XX,"assign") <- attr(X,"assign") ## possible bug here
            return(XX)
            
        }
    }else if(test.lme || test.gls){
        FUN.model.matrix <- function(model, formula, data){
            model.matrix(formula, data = data, fixed.only = FALSE)
        }
    }else if(test.lm){
        FUN.model.matrix <- function(model, formula, data){
            model.matrix(model, data = data, fixed.only = FALSE)
        }
    }

    FUN.formula <- function(model){
        ff <- try(formula(model), silent = TRUE)
        if("try-error" %in% class(ff)){
            ff <- lavaSearch2_evalInParentEnv(model$call[[2]])
        }
        if("formula" %in% class(ff) == FALSE){
            stop("Unable to extract the formula from the model \n")
        }
        return(ff)
    }
    
    if(missing(FUN.df)){
        if(test.lmer || test.gls){
            FUN.df <- function(model, level){
                NULL
            }
        }else if(test.lme){
            FUN.df <- function(model, level){
                min(model$fixDF$X)
            }
        }else if(test.lm || test.gam){
            FUN.df <- function(model, level){
                model$df.residual
            }

        }
    }
    
    ## ** normalize input
    formula.model <- FUN.formula(model)
    interval <- match.arg(interval, c("confidence","prediction"))
    design.df <- as.data.table(lavaSearch2::extractData(model, design.matrix = FALSE))[,.SD,.SDcols=all.vars(formula.model)]
    
    design.numeric <- sapply(design.df, is.numeric)
    design.factor <- sapply(design.df, is.factor)

    design.NonNumeric <- names(which(design.numeric==FALSE))
    if(length(design.NonNumeric)>0){
        xlevels <- sapply(design.NonNumeric, function(x){
            list(levels(as.factor(design.df[[x]])))
        })
    }

    all.vars.model <- all.vars(formula.model)
    name.Y <- all.vars.model[1]
    name.X.df <- all.vars.model[-1]
    if (any(var %in% name.X.df == FALSE) ) {
        stop("Unknown variable(s): \"",paste(var[var %in% name.X.df == FALSE], collapse = "\" \""),"\" \n",
             "argument \'var\' should be in \"",paste(name.X.df,collapse = "\" \""),"\" \n")
    }
    if(length(keep.intercept)==1){
        keep.intercept <- rep(keep.intercept,2)
    }
 
    ## ** Reference level
    n <- NROW(design.df)
    p <- length(name.X.df)
    newdata.fit <- as.data.frame(design.df)

    var.set0 <- names(newdata.fit)[which(names(newdata.fit) %in% var == TRUE)]
    if(length(var.set0)){
        for(i0 in var.set0){ # i0 <- var.set0[1]
            if(design.numeric[i0]){
                newdata.fit[[i0]][] <- 0
            }else{
                newdata.fit[[i0]] <- factor(xlevels[[i0]][1], levels = xlevels[[i0]])
            }
            
        }
    }
    
    ## ** intercept
    beta <- FUN.coef(model)
    design.mat <- FUN.model.matrix(model, formula = formula.model, data = newdata.fit)    
    test.intercept <- as.logical(attr(terms(formula.model),"intercept"))
    if(test.intercept){
        name.intercept <- "(Intercept)"
    }else if(any(design.factor)){ # find reference level
        attr.assign <- attr(design.mat, "assign")
        indexRef <- which(tapply(attr.assign, attr.assign, function(x){any(duplicated(x))}))
        name.intercept <- names(beta)[match(indexRef, attr.assign)]
    }else{
        name.intercept <- NULL
    }

    ## ** partial fit
    if(!is.null(name.intercept) && keep.intercept[1]){
        ## add the residual due to the intercept
        design.mat[,name.intercept] <- 0
    }

    newdata.index <- try(as.numeric(rownames(design.mat)),silent=FALSE)
    if(!is.null(newdata.index) && !inherits(newdata.index,"try-error")){
        index.nna <- newdata.index
    }else{
        index.nna <- setdiff(1:NROW(newdata.fit),model$na.action)
    }
    design.df$pFit <- as.numeric(NA)
    design.df$pFit[index.nna] <- as.numeric(design.mat %*% beta)

    ## ** random effect
    if(conditional){ # https://stackoverflow.com/questions/25538199/design-matrix-for-mlm-from-librarylme4-with-fixed-and-random-effects
        design.df$ranef[index.nna] <- FUN.ranef(model)
    }else{
        design.df$ranef <- 0
    }

    ## ** partial residuals
    design.df$pResiduals <- newdata.fit[[name.Y]] - design.df$pFit - design.df$ranef
    
    ## ** partial fitted values
    ## newdata
    ls.forGrid <- lapply(name.X.df, function(x){ ## x <- name.X.df[1]
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
            return(intersect(unique(design.df[[x]]),xlevels[[x]]))
        }
    })

    ## create the data.frame
    grid.predict <- expand.grid(ls.forGrid, KEEP.OUT.ATTRS = FALSE)
    names(grid.predict) <- name.X.df
    grid.predict[[name.Y]] <- 0
    ## restaure factor variables
    if(any(design.numeric == FALSE)){         
        for(iVar in names(which(design.numeric == FALSE))){
            grid.predict[[iVar]] <- factor(grid.predict[[iVar]], levels = xlevels[[iVar]])
        }
    }
    ## convert to design matrix
    design.grid.predict <- FUN.model.matrix(model, formula = formula.model, grid.predict)
    ## remove intercept
    if(!is.null(name.intercept) && keep.intercept[1] == FALSE){
        ## remove the intercept from the prediction
        ## when it is not among the variable of interest
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

    design.grid.predict2 <- design.grid.predict
    if(!is.null(name.intercept) && keep.intercept[2] == FALSE){
        ## remove the intercept from the prediction
        ## when it is not among the variable of interest
        design.grid.predict2[,name.intercept] <- 0
    }
    se.tempo <- apply(design.grid.predict2, 1, function(x){
        sqrt(rbind(x) %*% model.vcov %*% cbind(x) + model.sigma2)
    })
    ## sqrt(diag(model.vcov)[2])
    ## se.tempo <- summary(model)$coef[2,2]
    
    grid.predict$fit.lower <- grid.predict$fit - z * se.tempo
    grid.predict$fit.upper <- grid.predict$fit + z * se.tempo        

    out <- as.data.table(design.df)
    attr(out, "partialFit") <- as.data.table(grid.predict)
    attr(out, "var") <- var
    attr(out, "interval") <- interval
    attr(out, "level") <- level
    attr(out, "name.Y") <- name.Y
    class(out) <- append("partialResiduals",class(out))
    return(out)
}


## * predict functions
#' @title prediction function for partialResiduals
#' @description Suggestion of functions to obtain confidence intervals with mixed models
#' @name FUN.predict
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

## * predict_merTools
#' @rdname FUN.predict
#' @export
predict_merTools <- function(object, newdata, level,conditional, interval, ...){
    if(missing(interval)){interval <- FALSE}
    if(missing(level)){level <- 0.95}
    requireNamespace("merTools")
    res <- merTools::predictInterval(object, newdata = newdata, level = level,
                                     type = "linear.prediction", n.sims = 1e3,
                                     which = if(conditional){"full"}else{"fixed"},
                                     include.resid.var = (interval=="prediction"))
    return(list(fit = res))
}

## * predict_AICcmodavg 
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

