### mcrGam.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 19 2018 (10:08) 
## Version: 
## Last-Updated: jun 20 2018 (09:30) 
##           By: Brice Ozenne
##     Update #: 12
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * mcrGam (documentation)
#' @title Monotone cr spline
#' @description Essentially group the line of codes in the example section of \code{mgcv::pcls}.
#' @name mcrGam
#' 
#' @param object [gam] output of the \code{mgcv::gam} function.
#' @param increasing [logical] should the function should be constrained to be increasing.
#' Otherwise it will be constrained to be decreasing.
#' @param lower [numeric] This specifies the lower bound on the spline unless it is NA in which case no lower bound is imposed.
#' @param upper [numeric] This specifies the upper bound on the spline unless it is NA in which case no upper bound is imposed.
#' @param ... Not used. For compatibility with the generic function.
#' 
#' @examples
#' if(require(mgcv)){
#' 
#' ## Generate data from a monotonic truth.
#' set.seed(10)
#' x <- runif(100)*4-1;x <- sort(x);
#' f <- exp(4*x)/(1+exp(4*x)); y <- f+rnorm(100)*0.1; plot(x,y)
#' dat <- data.frame(x=x,y=y)
#' 
#' ## Show regular spline fit (and save fitted object)
#' 
#' f.ug <- gam(y~s(x,k=10,bs="cr")); lines(x,fitted(f.ug))
#' 
#' ## Create Design matrix, constraints etc. for monotonic spline....
#' sm <- smoothCon(s(x,k=10,bs="cr"),dat,knots=NULL)[[1]]
#' F <- mono.con(sm$xp);   # get constraints
#' G <- list(X=sm$X,C=matrix(0,0,0),sp=f.ug$sp,p=sm$xp,y=y,w=y*0+1)
#' G$Ain <- F$A;G$bin <- F$b;G$S <- sm$S;G$off <- 0
#' p <- pcls(G);  # fit spline (using s.p. from unconstrained fit)
#' fv<-Predict.matrix(sm,data.frame(x=x))%*%p
#' lines(x,fv,col=2)
#'
#' ## Proposed function
#' montone.f.ug <- mcrGam(f.ug)
#' monotone.pred <- predict(montone.f.ug)
#' 
#' if(require(testthat)){
#' expect_equal(montone.f.ug$mcr.p,p)
#' expect_equal(montone.f.ug$mcr.sm,sm)
#' expect_equal(monotone.pred,fv)
#' }
#'
#' }
#' @export
mcrGam <- function(object, ...){
    UseMethod("mcrGam") 
}

## * mcrGam (code)
#' @method mcrGam gam
#' @rdname mcrGam
#' @export
mcrGam.gam <- function(object, increasing = TRUE, lower = NA, upper = NA, ...){

    ## ** extract information
    ## methods(class = class(object)[1])
    
    ## data
    data <- model.frame(object)

    ## response var
    formula.object <- formula(object)
    response.var <- all.vars(update(formula.object, ".~0"))
    
    ## spline term
    term.object <- terms(formula.object,
                         specials = "s")

    index.spline <- attr(term.object,"specials")$s
    if(is.null(index.spline)){
        stop("Could not find any spline term \n")
    }
    if(length(index.spline)>1){
        stop("Only work for one spline term \n")
    }

    txt.spline <- rownames(attr(term.object,"factors"))[index.spline]
    smooth.spline <- eval(parse(text = txt.spline))
    ## names(smooth.spline)
    ## class(smooth.spline)
    
    if(class(smooth.spline)!="cr.smooth.spec"){
        stop("Only handle spline using \"cr\" basis \n",
             "Specify s(, bs = CR) \n")
    }
    ## keep.term <-  c("term", "bs.dim", "fixed", "dim", "p.order", "by", "label", "xt", "id", "sp")
    ## smooth.spline <- object$smooth[[1]][keep.term]
    ## class(smooth.spline) <- "cr.smooth.spec"
    ## lapply(keep.term, function(index){
    ## identical(smooth.spline[[index]], eval(parse(text = txt.spline))[[index]])
    ## })

    ## ** Extract parametrisation of the splines
    sm <- mgcv::smoothCon(smooth.spline, data = data)[[1]]
    ## sm$xp knot location
    
    ## ** Find linear constrains sufficient for monotonicity
    mc <- mgcv::mono.con(sm$xp,
                         up = increasing,
                         lower = lower,
                         upper = upper) # monotonicity constraints
    
    M <- list(X = sm$X, y = data[[response.var]], #design matrix, outcome
              C = matrix(0,0,0), #equality constraints (none)
              Ain = mc$A, bin = mc$b, #inequality constraints
              sp = object$sp, p = sm$xp, #initial guesses for param estimates
              S = sm$S, #smoothness penalty matrix
              w = data[[response.var]]*0+1, off=0 #weights, offset
              )

    ## ** fit spine using penalized constrained least squares
    p <- mgcv::pcls(M)

    ## ** export
    object$mcr.sm <- sm
    object$mcr.p <- p
    class(object) <- "mcrGam"
    return(object)
}

## * print.mcrGam
#' @method print mcrGam
#' @export
print.mcrGam <- function(x, ...){
    cat("Spline coefficients: \n")
    print(x$mcr.p)
}

## * predict.mcrGam
#' @method predict mcrGam
#' @export
predict.mcrGam <- function(object, newdata = NULL, ...){

    ## extract data
    if(is.null(newdata)){
        data <- model.frame(object)
    }
    data.table::setDT(data)

    ## find spline variable
    formula.object <- formula(object)
    term.object <- terms(formula.object,
                         specials = "s")
    index.spline <- attr(term.object,"specials")$s
    spline.var <- all.vars(formula.object)[index.spline]

    ## ** Predict values 
    return(mgcv::Predict.matrix(object =object$mcr.sm,
                                data = data[,.SD,.SDcols = spline.var]) %*% object$mcr.p)
}



######################################################################
### mcrGam.R ends here
