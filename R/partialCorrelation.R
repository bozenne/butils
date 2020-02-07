### partialCorrelation.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb  6 2020 (13:28) 
## Version: 
## Last-Updated: feb  6 2020 (17:37) 
##           By: Brice Ozenne
##     Update #: 82
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
    function(object, var, fisher.transform, cluster) UseMethod("partialCorrelation")

## * partialCorrelation.lm
##' @rdname partialCorrelation
##' @export
partialCorrelation.lm <- function(object, var, fisher.transform = FALSE, cluster = NULL){

    ## ** normalize arguments
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

    ## ** remove covariate effect
    if(length(X)>0){
        object1.formula <- paste0(Y,"~",paste(X,collapse="+"))
        object2.formula <- paste0(var,"~",paste(X,collapse="+"))
    }else{
        object1.formula <- paste0(Y,"~1")
        object2.formula <- paste0(var,"~1")
    }
    object1 <- update(object, formula = as.formula(object1.formula))
    object2 <- update(object, formula = as.formula(object2.formula))

    ## ** compute residuals and their correlation
    n.obs <- nobs(object) + length(object$na.action)
    df.res <- data.frame(matrix(NA, nrow = n.obs, ncol = 2,
                                dimnames = list(NULL, c("res1","res2"))))
    df.res[setdiff(1:n.obs,object1$na.action),"res1"] <- scale(residuals(object1))
    df.res[setdiff(1:n.obs,object2$na.action),"res2"] <- scale(residuals(object2))
    if(!is.null(cluster)){
        df.res[[cluster]] <- NA
        df.res[setdiff(1:n.obs,object1$na.action),cluster] <- extractData(object1)[[cluster]]
        df.res[setdiff(1:n.obs,object2$na.action),cluster] <- extractData(object2)[[cluster]]
    }
    e.lmres <- lm(res1 ~ res2, data = df.res)

    ## ** extra info for export
    e.lmres <- sCorrect(e.lmres)
    e.lmres$sCorrect$data
    out <- summary2(e.lmres)$table2["res2",]
    rownames(out) <- NULL
    e.lmres.iid <- iid2(e.lmres, robust = FALSE, cluster = cluster)
    attr(out,"iid") <- setNames(e.lmres.iid[,"res2"],attr(e.lmres.iid,"cluster"))
    ## sqrt(sum(attr(out,"iid")^2))

    if(fisher.transform){ ## to match lava::partialcor
        out.trans <- transformSummaryTable(out, transform = "atanh")
        out.trans$std.error <- sqrt(1/(stats::nobs(e.lmres) - length(X) - 1 - 3))
        out[,"ci.lower"] <- tanh(out.trans$estimate + qnorm(0.025) * out.trans$std.error) ##tanh(out.trans[,c("ci.lower","ci.upper")])
        out[,"ci.upper"] <- tanh(out.trans$estimate + qnorm(0.975) * out.trans$std.error) ##tanh(out.trans[,c("ci.lower","ci.upper")])
        out[,"df"] <- NA
        out[,"p.value"] <- 2*(1-pnorm(abs(out.trans$estimate)/out.trans$std.error))
    }
    return(out)
}

## * partialCorrelation.mmm
partialCorrelation.mmm <- function(object, var, fisher.transform = FALSE, cluster = NULL){

    ## ** compute partial correlation for each model
    ls.out <- lapply(object, function(iM){ ## iM <- object[[5]]
        iVar <- grep(var, all.vars(formula(iM)), value = TRUE)
        if(length(iVar)>0){
           return(partialCorrelation(iM, var = iVar, fisher.transform = fisher.transform, cluster = cluster))
        }else{
            return(NULL)
        }
    })
    test.null <- sapply(ls.out, length)
    if(all(test.null==0)){return(NULL)}
    ls.out <- ls.out[test.null>0]
    dt.out <- do.call(rbind,ls.out)

    ## ** collect iid across clusters
    n.model <- length(ls.out)
    name.model <- names(ls.out)
    if(!is.null(cluster)){
        unique.cluster <- unique(unlist(lapply(ls.out, function(iM){
            names(attr(iM,"iid"))
        })))
        n.cluster <- length(unique.cluster)
        iid.out <- matrix(NA, nrow = n.cluster, ncol = n.model ,
                                  dimnames = list(unique.cluster, name.model))
        for(iM in 1:n.model){ ## iM <- 1
            iid.out[names(attr(ls.out[[iM]],"iid")),iM] <- as.double(attr(ls.out[[iM]],"iid"))
        }
    }else{
        iid.out <- do.call(cbind,lapply(ls.out, function(iM){
            attr(iM,"iid")
        }))
    }

    ## ** compute global vcov
    vec.sd <- apply(iid.out,2,function(iCol){sqrt(sum(iCol^2,na.rm = TRUE))})
    M.cor <- cor(iid.out, use = "pairwise")
    vcov.out <- M.cor * tcrossprod(vec.sd)
    ## sqrt(diag(attr(out,"vcov")))-out$std.error

    ## ** convert to multcomp
    linfct <- matrix(0, nrow = n.model, n.model,
                     dimnames = list(name.model, name.model))
    diag(linfct) <- 1
    out <- list(model = NULL,
                linfct = linfct,
                rhs = rep(0, n.model),
                coef = dt.out$estimate,
                vcov = vcov.out,
                df = if(any(!is.na(dt.out$df))){dt.out$df}else{0},
                alternative = "two.sided",
                type = NULL,
                robust = FALSE)
    class(out) <- c("glht2","glht")
    return(out)
}

######################################################################
### partialCorrelation.R ends here
