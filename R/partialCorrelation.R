### partialCorrelation.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb  6 2020 (13:28) 
## Version: 
## Last-Updated: feb  6 2020 (14:27) 
##           By: Brice Ozenne
##     Update #: 24
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * partialCorrelation
##' @title Compute the Partial Correlation Between the Outcome and a Covariate
##' @description Compute the partial correlation between the outcome and a covariate in a linear model.
##' @name partialCorrelation
##'
##' @param object an \code{lm} object.
##' @param var [character] the covariate with which the partial correlation should be computed.
##' @param fisher.transform [logical] should the p-value/confidence intervals be computed using Fisher's Z transform.
##' Otherwise a Student's t-distribution is used
##'
##' @examples
##' n <- 1e2
##' set.seed(10)
##' df <- data.frame(Y = rnorm(n),
##'                  X = rnorm(n),
##'                  K = as.character(rbinom(n, size = 3, prob = 0.5)))
##'
##' ## 1 covariate
##' e1.lm <- lm(Y~X, data = df)
##' partialCorrelation(e1.lm, var = "X")
##' cor.test(df$Y,df$X) ## same p-value, different CI
##'
##' ## 2 covariates
##' e2.lm <- lm(Y~X+K, data = df)
##' partialCorrelation(e2.lm, var = "X")
##' 
##' @export
`partialCorrelation` <-
    function(object, var, fisher.transform) UseMethod("partialCorrelation")

## * partialCorrelation.lm
##' @rdname partialCorrelation
##' @export
partialCorrelation.lm <- function(object, var, fisher.transform = FALSE){
    if(length(var)!=1){
        stop("Argument \'var\' must have length 1 \n")
    }
    object.formula <- formula(object)
    if(attr(terms(object.formula),"response")!=1){
        stop("Note implemented for \"lm\" object with multiple response variables \n")
    }
    name.var <- all.vars(object.formula)
    Y <- name.var[1]
    if(any(var %in% name.var == FALSE)){
        stop("Argument \'var\' must correspond to variables in the model \n")
    }

    if(var %in% names(object$xlevels)){
        stop("Argument \'var\' must correspond to a numeric variable \n")
    }
    X <- setdiff(name.var, c(Y,var))

    if(length(X)>0){
        object2.formula <- paste0(var,"~",Y,"+",paste(X,collapse="+"))
    }else{
        object2.formula <- paste0(var,"~",Y)
    }
    object2 <- update(object, formula = object2.formula)

    df.res <- data.frame(res1 = scale(residuals(object)),
                         res2 = scale(residuals(object2)))
    e.lmres <- lm(res1 ~ res2, data = df.res)
    e.lmres <- sCorrect(e.lmres)
    out <- summary2(e.lmres)$table2["res2",]
    attr(out,"iid") <- iid2(e.lmres)[,"res2"]

    ## cor.test(df.res$res1,df.res$res2)
    ## unlist(riskRegression::transformCIBP(estimate = out$estimate, se = cbind(out$std.error), null = 0, type = "atanh", band = FALSE, ci = TRUE, conf.level = 0.95,
    ##                               min.value = -1, max.value = 1, p.value = TRUE))
    if(fisher.transform){
        out.trans <- transformSummaryTable(out, transform = "atanh")
        out[,"ci.lower"] <- tanh(out.trans$estimate + qnorm(0.025) * out.trans$std.error) ##tanh(out.trans[,c("ci.lower","ci.upper")])
        out[,"ci.lower"] <- tanh(out.trans$estimate + qnorm(0.025) * out.trans$std.error) ##tanh(out.trans[,c("ci.lower","ci.upper")])
        out[,"df"] <- NA
        out[,"p.value"] <- 2*(1-pnorm(abs(out.trans$estimate)/out.trans$std.error))
    }
    return(out)
}

######################################################################
### partialCorrelation.R ends here
