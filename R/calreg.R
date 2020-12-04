### calreg.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 14 2020 (17:23) 
## Version: 
## Last-Updated: feb 20 2020 (17:01) 
##           By: Brice Ozenne
##     Update #: 156
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


## * calreg (documentation)
#' @title Regression On Unobserved Exposure Using Calibration.
#' @description Perform a linear regression on an unobserved exposure (X) using a proxy (Z) whose relationship with the exposure has been studied using an external dataset.
#'
#' @param formula formula for the linear model.
#' @param data [data.frame] dataset used to fit the linear model relating the outcome (Y) to the (unobserved) exposure (X).
#' @param calibration a \code{lm} object or \code{nls} object relating a proxy (Z) to the unobserved exposure (X).
#' @param method [character] Can be \code{"delta"} or \code{"MI"} to use, respectively, a delta method or multiple imputation.
#' @param n.impute [integer, >0] Number of imputed dataset to be used. Only relevant when \code{method="MI"}.
#'
#' @details Consider a first sample \eqn{(X_i,Z_i)_{i \in \{1,\ldots,m\}}} that is used to estimate \eqn{\alpha} in:
#' \deqn{X = f(alpha,Z) + \varepsilon_{\alpha}}
#' This is the model to give to the argument \code{calibration}.
#'
#' The aim is to use a second sample \eqn{(Y_j,Z_j)_{j \in \{1,\ldots,n\}}} to estimate \eqn{\beta_1} in:
#' \deqn{Y = \beta_0 + \beta_1 X + \varepsilon_{\beta}}
#' The formula of this model should be given to the argument \code{formula} and the dataset to the argument \code{data}.
#' 
#' The exposure \(X\) in the second sample is computed:
#'\itemize{
#' \item based on the conditional expectation of the exposure given the proxy from the first model (\code{method="delta"}).
#' \item based on multiple sampling of the coefficients from the first model (\code{method="MI"}).
#' For each sample an exposure is computed, a linear model is then estimated based on this exposure. The results are then pooled using \code{mice::pool}.
#' }
#' 
#' When using the delta method, the uncertainty is decomposed into two parts:
#' \itemize{
#' \item one related to the finite number of observations in the second sample.
#' \item one related to the estimation of the parameters in the calibration model, to account for the fact that \(X\) is estimated and not observed.
#' }
#' 
#' @return A data.frame containing the estimates, standard errors, confidence intervals and p-values for each regression coefficient.
#' The output has an attribute \code{"regression"} containing the fitted linear model (ignoring the uncertainty related to the calibration)
#' and an attribute \code{"var.add"} representing additional variance-covariance matrix due to the calibration.
#'
#' @examples
#' library(lava)
#' n <- 1e2
#' 
#' ## linear case
#' mSim.lin <- lvm(fMRI ~ occ, occ ~ blood)
#' distribution(mSim.lin, ~blood) <- uniform.lvm(-0.9,2)
#'
#' set.seed(10)
#' d1.lin <- sim(mSim.lin, n = n)[,c("occ","blood"),drop=FALSE]
#' d2.lin <- sim(mSim.lin, n = n)[,c("fMRI","blood"),drop=FALSE]
#'
#' e1.lin <- lm(occ ~ blood, data = d1.lin)
#' res.lin <- calreg(fMRI ~ occ, data = d2.lin, calibration = e1.lin)
#' res.lin
#' summary(attr(res.lin, "regression"))$coef
#'
#' ## non-linear case
#' mSim.nlin <- lvm(fMRI ~ occ, occ[mu:0.1] ~ 0*blood)
#' distribution(mSim.nlin, ~blood) <- uniform.lvm(-0.9,3)
#' constrain(mSim.nlin, mu~blood) <- function(x){2*x/(1+x)}
#'
#' set.seed(10)
#' d1.nlin <- sim(mSim.nlin, n = n)[,c("occ","blood"),drop=FALSE]
#' d2.nlin <- sim(mSim.nlin, n = n)[,c("fMRI","blood"),drop=FALSE]
#'
#' ## gg <- ggplot(d1.nlin, aes(x = blood)) + geom_point(aes(y = occ))
#'
#' e1.nlin <- nls(occ ~ (occmax * blood)/(EC + blood), data = d1.nlin,
#'          start = list(occmax = 1, EC = 1))
#' d1.nlin$fit <- fitted(e1.nlin)
#' 
#' ## gg + geom_line(data = d1.nlin, aes(y = fit), color = "red")
#'
#' res.nlin <- calreg(fMRI ~ occ, data = d2.nlin, calibration = e1.nlin)
#' res.nlin
#' summary(attr(res.nlin, "regression"))$coef
#' 

## * calreg (code)
##' @rdname calreg
##' @export
calreg <- function(formula, data, fitter = "lm", calibration, method = "delta", n.impute = 50){

    ## ** extract information
    formula.calibration <- formula(calibration)
    name.X <- all.vars(formula.calibration)[1]
    args.txt <- names(coef(calibration))

    alpha <- coef(calibration)
    name.alpha <- names(alpha)
    vcov_alpha <- vcov(calibration)

    
    ## ** check arguments
    if(!inherits(calibration,"lm") && !inherits(calibration,"nls")){
        stop("Only compatible with lm and nls objects \n")
    }
    method <- match.arg(method, choices = c("delta","MI"))

    if(name.X %in% all.vars(formula)[-1] == FALSE){
        stop("The outcome of the calibration model should be an covariate in the argument \'formula\' \n")
    }

    fitter <- match.arg(fitter, c("lm","lmer"))
    getCoef <- switch(fitter,
                      "lm" = coef,
                      "lmer" = lme4::fixef)

    ## ** prepare prediction
    if(inherits(calibration,"lm")){
        Z <- model.matrix(calibration)
        refit <- function(p){
            if(is.null(p)){
                p <- coef(calibration)
            }
            data[[name.X]] <- Z %*% p
            return(do.call(fitter, args = list(formula = formula, data = data)))
        }
    }else if(inherits(calibration,"nls")){
        class(calibration) <- append("nls.calreg",class(calibration))
        refit <- function(p){
            data[[name.X]] <- predict(calibration, newdata = data, newparam = p)
            return(do.call(fitter, args = list(formula = formula, data = data)))
        }
    }

    ## ** fit model based on the expected exposure conditional on the proxy
    regression <- refit(NULL) ## default values (same as refit(alpha))
    Sregression <- summary(regression)

    if(any(is.na(getCoef(regression)))){
        warning("One of the regression parameter could not be estimated.\n")
        table2 <- data.frame(matrix(NA, nrow = nrow(Sregression$coef), ncol = 7,
                                    dimnames = list(rownames(Sregression$coef), c("estimate","std.error","df","ci.lower","ci.upper","statistic","p.value"))))
        attr(table2,"regression") <- regression
        attr(table2,"add.var") <- NA
        return(table2)
    }

    ## ** approach 1: delta method
    if(method == "delta"){
    
        ## apply delta method
        ## For f: Rn -> Rm calculate the m,n Jacobian dy/dx (i.e. numerical approxiamtion of the first derivative)
        dcoef <- numDeriv::jacobian(func = function(x){getCoef(refit(x))},
                                    x = alpha)
        var.add <- dcoef %*% vcov_alpha %*% t(dcoef)
        
        ## dcoef[1,,drop=FALSE] %*% vcov_alpha %*% t(dcoef[1,,drop=FALSE])
        ## dcoef[2,,drop=FALSE] %*% vcov_alpha %*% t(dcoef[2,,drop=FALSE])
        ## numDeriv::grad(func = function(x){getCoef(refit(x))}[1],
        ##                    x = alpha)
        ## numDeriv::grad(func = function(x){getCoef(refit(x))}[2],
        ##                    x = alpha)
        
        ## dcoef[1,,drop=FALSE] %*% vcov_alpha %*% t(dcoef[1,,drop=FALSE]) - var.add[1,1]
        newcov <- as.matrix(vcov(regression)) + var.add ## as.matrix necessary for lmer object to convert sparse matrix to matrix
        table2 <- data.frame(matrix(NA, nrow = nrow(Sregression$coef), ncol = 7,
                                    dimnames = list(rownames(Sregression$coef), c("estimate","std.error","df","ci.lower","ci.upper","statistic","p.value"))))
        table2$estimate <- as.double(Sregression$coef[,"Estimate"])
        table2$std.error <- as.double(sqrt(diag(newcov)))
        if(inherits(regression,"lmerModLmerTest")){
            table2$df <- unname(summary(regression)$coef[,"df"])
        }else if(inherits(regression,"lm")){
            table2$df <- regression$df.residual
        }else{
            table2$df <- Inf
        }
        table2$ci.lower <- table2$estimate + qt(0.025, df = table2$df) * table2$std.error
        table2$ci.upper <- table2$estimate + qt(0.975, df = table2$df) * table2$std.error
        table2$statistic <- table2$estimate/table2$std.error
        table2$p.value <- 2*(1-pt(abs(table2$statistic), df = table2$df))

        attr(table2, "regression") <- regression
        attr(table2, "var.add") <- var.add
    }

    ## ** approach 2: multiple imputation
    if(method %in% "MI"){
        M.impute <- MASS::mvrnorm(n = n.impute, mu = alpha, Sigma = vcov_alpha)
        lmer <- lme4::lmer
        regression.MI <- apply(M.impute, 1, function(iCoef){
            refit(iCoef)
        })
        names(regression.MI) <- paste0("M",1:length(regression.MI))
        pool.MI <- mice::pool(mice::as.mira(regression.MI))
        table2 <- summary(pool.MI)
        table2$ci.lower <- table2$estimate + qt(0.025, df = table2$df) * table2$std.error
        table2$ci.upper <- table2$estimate + qt(0.975, df = table2$df) * table2$std.error
        table2 <- table2[,c("estimate","std.error","df","ci.lower","ci.upper","statistic","p.value")]

        attr(table2, "regression") <- regression
        attr(table2, "var.add") <- table2[name.X,"std.error"]^2 - vcov(regression)[name.X,name.X]
    }

    ## ** export
    return(table2)
    
}

## * predict.nls.calreg
## prediction with new coefficient value 
predict.nls.calreg <- function(object, newdata = NULL, newparam = NULL, se.fit = FALSE, ...){
    class(object) <- setdiff(class(object),"nls.calreg")

    if(!is.null(newparam)){
        if(se.fit){
            stop("Argument \'se.fit\' must be FALSE when argument \'newparam\' is not NULL. \n")
        }

        ## ** check names
        args.param <- names(object$m$getPars())
        object.env <- object$m$getEnv()
        object.env.save <- mget(ls(envir = object.env), envir = object.env)
        if(is.null(names(newparam)) || (any(sort(names(newparam)) != sort(args.param)))){
            stop("Incorrect specification of argument \'newparam\' \n",
                 "Must be a vector/list of numeric values with name(s): \"",paste(args.param, collapse = "\" \""),"\"\n")
        }

        if(any(args.param %in% names(object.env) == FALSE)){ ## then parameters are specified via a vector
            args.param <- unique(gsub("\\d+$", replacement = "", x = args.param))
            if(length(args.param)!=1){
                stop("Could not identify a unique vector of parameters \n")
            }
            newparam <- setNames(list(as.double(newparam)), args.param)
        }else{            
            newparam <- as.list(newparam)
        }
        
        ## ** assign param
        n.param <- length(args.param)
        for(iP in 1:n.param){
            iParam <- args.param[iP]
            base::assign(x = iParam, value = as.numeric(newparam[[iParam]]), envir = object.env)
        }
    }

    out <- predict(object, newdata = newdata, se.fit = se.fit, ...)

    if(!is.null(newparam)){
        ## ** restore param
        for(iP in 1:n.param){
            iParam <- args.param[iP]
            base::assign(x = iParam, value = object.env.save[[iParam]], envir = object.env)
        }
    }
        
    return(out)
}

######################################################################
### calreg.R ends here
