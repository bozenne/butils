### officer.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 18 2019 (09:35) 
## Version: 
## Last-Updated: jan 23 2020 (18:30) 
##           By: Brice Ozenne
##     Update #: 51
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * multcomp2table - documentation
##' @title Gather P-values and Confidence Intervals in a Table.
##' @description Gather p-values and confidence intervals in a table.
##' @name multcomp2table
##'
##' @param object a fitted model.
##' @param link [character vector] null hypotheses or coefficients to be tested.
##' @param transform [function] function to backtransform the estimates and the associated confidence intervals.
##' @param seed [integer] if not \code{NULL}, integer used to initialize the random number generator. See \code{\link{set.seed}} for details.
##' @param method.glht [character or function] function used to extract the coefficients and variance-covariance matrix from the object.
##' Recommanded: \code{"glht"} from the multcomp package or \code{"glht2"} from the lavaSearch2 package.
##' @param conf.level [numeric, 0-1] Confidence level of the interval.
##' @param method.multcomp [character] the method used to adjust the p-value and confidence intervals (CIs) for multiplicity.
##' Note that simultaneous CIs are available only for \code{method.multcomp="none"},  \code{method.multcomp="bonferroni"}, and  \code{method.multcomp"single-step"}.
##' @param digits [integer] if not \code{NULL}, the number of decimal places used for the estimate and the associated confidence intervals.
##' @param digits.p.value [integer] if not \code{NULL}, the number of decimal places used for the p-values.
##' @param ... arguments passed to \code{method.glht}.
##' 
##' 

## * multcomp2table - examples
##' @rdname multcomp2table
##' @examples
##' if(require(multcomp)){
##'
##' m <- lvm(Y~X1+X2)
##' d <- sim(m, n = 100)
##'
##' ## lm object
##' e.lm <- lm(Y~X1+X2, data = d)
##' multcomp2table(e.lm, link = c("X1=0","X2=1"))
##' multcomp2table(e.lm, link = c("X1=0"))
##'
##' ## gls object
##' if(require(nlme)){
##' e.gls <- gls(Y~X1+X2, data = d)
##' multcomp2table(e.gls, link = c("X1=0","X2=1"))
##' multcomp2table(e.gls, link = c("X1=0"))
##' }
##' 
##' ## lvm object
##' if(require(lava) & require(lavaSearch2)){
##' e.lvm <- estimate(m, data = d)
##' multcomp2table(e.lvm, method.glht = "glht2", link = c("Y~X1"))
##' multcomp2table(e.lvm, method.glht = "glht2", link = c("Y~X1","Y~X2"))
##' multcomp2table(e.lvm, method.glht = "glht2", link = c("Y~X1","Y~X2"), rhs = c(1,1))
##' }
##' 
##' }

## * multcomp2table - code
##' @rdname multcomp2table
##' @export
multcomp2table <- function(object, link, transform = function(x){x}, seed = NULL,
                           method.glht = "glht", conf.level = 0.95, method.multcomp = NULL,
                           digits = 3, digits.p.value = 3, ...){
    if(!is.null(seed)){
        set.seed(NULL)
    }
    
    ## ** normalization
    alpha <- 1-conf.level
    if(inherits(object,"glht")){
        if(missing(link)){
            link <- rownames(object$linfct)
            e.glht <- object
        }else {
            stop("Argument \'link\' should not be specified when using a \"glht\" object \n")
        }
    }else{
        ## check object
        object.coef <- names(coef(object))
        vec.test <- sapply(link, grepl, pattern = paste0(object.coef,collapse="|"))
        if(any(vec.test == FALSE)){
            txt <- link[vec.test == FALSE]
            stop("Incorrect specification of argument \'link\' \n",
                 "links no found in object: \"",paste(txt, collpase = "\" \""),"\"\n")
        }
        e.glht <- do.call(method.glht, args = list(object, linfct = link, ...))
    }
    if(!is.null(method.multcomp)){
        eS.glht <- summary(e.glht, test = multcomp::adjusted(method.multcomp))
    }else if(!is.null(e.glht$test)){
        eS.glht <- e.glht
        method.multcomp <- e.glht$test$type
    }else{
        stop("The argument \'method.multcomp\' need to be specified \n",
             "Consider, for example, setting it to \"single-step\" \n")
    }
    
    
    ## ** prepare output
    n.link <- length(link)
    out <- as.data.frame(matrix(NA, nrow = n.link, ncol = 6,
                                dimnames = list(1:n.link,c("Parameter","Estimate","Std. Err.","Lower","Upper","p.value"))))

    out[,"Parameter"] <- link

    ## ** extract p-value and standard errors
    out[,c("Std. Err.","p.value")] <- do.call(cbind, eS.glht$test[c("sigma","pvalues")])

    ## ** extract estimates and CIs
    if(method.multcomp == "single-step"){
        out[,c("Estimate","Lower","Upper")] <- stats::confint(eS.glht, level = conf.level, calpha = adjusted_calpha())$confint
    }else if(method.multcomp == "none"){
        out[,c("Estimate","Lower","Upper")] <- stats::confint(eS.glht, level = conf.level, calpha = univariate_calpha())$confint
    }else if(method.multcomp == "bonferroni"){
        out[,c("Estimate","Lower","Upper")] <- stats::confint(eS.glht, level = 1-(1-conf.level)/NROW(eS.glht$linfct), calpha = univariate_calpha())$confint
    }else{
        out[,"Estimate"] <- eS.glht$test$coefficients
        out[,"Lower"] <- as.numeric(out[,"Lower"])
        out[,"Upper"] <- as.numeric(out[,"Upper"])
    }

    ## ** back transform
    out[,"Estimate"] <- transform(out[,"Estimate"])
    out[,"Upper"] <- transform(out[,"Upper"])
    out[,"Lower"] <- transform(out[,"Lower"])
    
    ## ** round
    if(!is.null(digits)){
        out[,c("Estimate","Std. Err.","Lower","Upper")] <- round(out[,c("Estimate","Std. Err.","Lower","Upper")], digits = digits)
        out[,"p.value"] <- format.pval(out["p.value"], digits =  digits.p.value, eps = 10^{-digits.p.value})
    }

    ## ** export
    return(out)
}


######################################################################
### officer.R ends here
