### officer.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 18 2019 (09:35) 
## Version: 
## Last-Updated: feb 18 2019 (10:21) 
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

## * multcomp2table - documentation
##' @title Gather P-values and Confidence Intervals in a Table.
##' @description Gather p-values and confidence intervals in a table.
##' @name multcomp2table
##'
##' @param object a fitted model.
##' @param link [character vector] null hypotheses or coefficients to be tested.
##' @param method.multcomp [character or function] function used to extract the coefficients and variance-covariance matrix from the object.
##' Recommanded: \code{"glht"} from the multcomp package or \code{"glht2"} from the lavaSearch2 package.
##' @param conf.level [numeric, 0-1] Confidence level of the interval.
##' @param adjust.multiple.comparison [logical] Should p-values and confidence intervals be adjusted for multiple testing using single step Dunnett.
##' @param digits [integer] if not \code{NULL}, the number of decimal places.
##' @param ... arguments passed to \code{method.multcomp}.
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
##' multcomp2table(e.lvm, method.multcomp = "glht2", link = c("Y~X1"))
##' multcomp2table(e.lvm, method.multcomp = "glht2", link = c("Y~X1","Y~X2"))
##' multcomp2table(e.lvm, method.multcomp = "glht2", link = c("Y~X1","Y~X2"), rhs = c(1,1))
##' }
##' 
##' }

## * multcomp2table - code
##' @rdname multcomp2table
##' @export
multcomp2table <- function(object, link,
                           method.multcomp = "glht", conf.level = 0.95, adjust.multiple.comparison = TRUE,
                           digits = 3, ...){

    ## check
    object.coef <- names(coef(object))
    vec.test <- sapply(link, grepl, pattern = paste0(object.coef,collapse="|"))
    if(any(vec.test == FALSE)){
        txt <- link[vec.test == FALSE]
        stop("Incorrect specification of argument \'link\' \n",
             "links no found in object: \"",paste(txt, collpase = "\" \""),"\"\n")
    }

    alpha <- 1-conf.level
    
    ## run
    n.link <- length(link)
    out <- as.data.frame(matrix(NA, nrow = n.link, ncol = 6,
                                dimnames = list(1:n.link,c("Parameter","Estimate","Std. Err.","Lower","Upper","p.value"))))

    e.glht <- do.call(method.multcomp, args = list(object, linfct = link, ...))
    out[,"Parameter"] <- link

    if(adjust.multiple.comparison){
        out[,c("Estimate","Lower","Upper")] <- stats::confint(e.glht, level = conf.level)$confint
        out[,c("Std. Err.","p.value")] <- do.call(cbind, summary(e.glht, test = multcomp::adjusted("single-step"))$test[c("sigma","pvalues")])
    }else{
        out[,c("Estimate","Lower","Upper")] <- stats::confint(e.glht, calpha = qt(1-alpha/2, df = e.glht$df))$confint
        out[,c("Std. Err.","p.value")] <- do.call(cbind, summary(e.glht, test = multcomp::univariate())$test[c("sigma","pvalues")])
    }
    
    if(!is.null(digits)){
        out[,c("Estimate","Std. Err.","Lower","Upper")] <- round(out[,c("Estimate","Std. Err.","Lower","Upper")], digits = digits)
        out[,"p.value"] <- format.pval(out["p.value"], digits =  digits, eps = 10^{-digits})
    }
    return(out)
}


######################################################################
### officer.R ends here
