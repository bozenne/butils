### R2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: aug  8 2017 (14:03) 
## Version: 
## last-updated: jan 27 2020 (15:20) 
##           By: Brice Ozenne
##     Update #: 106
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * calcR2 - documentation
#' @title Compute the R square using the predicted values
#' @description Compute the R square using the predicted values
#'
#' @name calcR2
#' 
#' @param model the model from which the R squared should be computed.
#' @param data the data that have been used to fit the model.
#' @param trace should the execution of the function be traced.
#' 
#' @details Compute the Generalized R2 defined by Buse (1973) using the formula:
#' \deqn{R^2 = 1- \frac{e S^{-1} e}{e0 S^{-1} e0}}
#' where iS is the inverse of the variance-covariance matrix of the observations.
#'       e is the vector of residuals of the full model.
#'       e0 is the vector of residuals of the reduced model.
#'       
#' Compute the Generalized R2 based on the log-likelihood:
#' \deqn{R^2 = 1 - exp(- LR / n)}
#' where n is the number of observations.
#'       LR is twice the log of the ratio of the likelihood between the full and reduced model.
#' 
#' Denoting X the covariate of interest, Y the outcome, and Z the other covariates, 
#' the partial R^2 is the same as:
#' \itemize{
#' \item the square of the partial correlation coefficient between X and Y, controlling for Z.
#' \item the partial \eqn{\eta^2} which is the ratio between the sum of squares from X (SSX) and the sum of SSX and the residual sum of squares (SSE)
#' }
#' 
#' @references
#' A. Buse (1973). Goodness of Fit in Generalized Least Squares Estimation. The American Statistician,, Vol. 27, No. 3.
#' Lonnie Magee (1990). R2 Measures Based on Wald and Likelihood Ratio Joint Significance Tests, The American Statistician, 44:3, 250-253.
#' https://stats.stackexchange.com/questions/64010/importance-of-predictors-in-multiple-regression-partial-r2-vs-standardized
#' 
#' @examples
#' library(lava)
#' library(nlme)
#' 
#' m <- lvm(Y~X1+X2+X3, G~1)
#' categorical(m, K=5, label = c("A","B","C","D","E")) <- ~G
#' d <- lava::sim(m, 1e2)
#' d <- d[order(d$G),]
#' d$Y <- d$Y + 0.5*as.numeric(as.factor(d$G))
#' 
#' ## linear model
#' m.lm <- lm(Y~X1+X2+X3, data = d)
#' calcR2(m.lm)
#' summary(m.lm)
#' 
#' summary(m.lm)$r.squared
#'
#' ff <- Y~X1+X2+X3
#' m.lm <- lm(ff, data = d)
#' calcR2(m.lm)
#' summary(m.lm)
#' 
#' dt <- as.data.table(d)
#' dt[,interaction := X1*X2]
#' m.lm2 <- lm(Y~X1+X2+X3+interaction, data = dt[G %in% c("A","B")])
#' calcR2(m.lm2)
#' 
#' ## gls model
#' m.gls1 <- gls(Y~X1+X2+X3, data = d)
#' calcR2(m.gls1)
#'
#' ## heteroschedasticity
#' m.gls2 <- gls(Y~X1+X2+X3,
#'               weights = varIdent(form = ~ 1 |G), data = d)
#' calcR2(m.gls2)
#'
#' ## correlation
#' m.gls2 <- gls(Y~X1+X2+X3,
#'               correlation = corCompSymm(form = ~ 1 |G), data = d)
#' calcR2(m.gls2)
#'

## * function - calcR2
#' @rdname calcR2
#' @export
calcR2 <- function(model, data = NULL, trace = FALSE){

    if(class(model) %in% c("lm","gls","lme") == FALSE){
        stop("Function only compatible with \'lm\', \'gls\', and \'lme\' objects \n")
    }
    
    ff <- formula(model)
    if(is.null(data)){
      data <- as.data.table(lavaSearch2::extractData(model, design.matrix = FALSE))
    }else{
      data <- copy(data.table::as.data.table(data))
    }
    n <- NROW(data)
    
    name.endo <- all.vars(ff[[2]])
    name.exo <- all.vars(ff[[3]])
    n.exo <- length(name.exo)
    
    ## V
    if(class(model)=="lm"){        
        iV <- NULL        
    }else if(class(model) %in% c("gls","lme")){
        std.residuals <-  attr(model$residuals,"std")
        if (is.null(model$modelStruct$corStruct)) { # no correlation matrix b
            V <- Matrix::Diagonal(x = std.residuals)
        }else{
            ls.V <- lapply(levels(model$groups), function(x){
                if(!inherits(model,"sCorrect")){
                    model <- sCorrect(model, ssc = model$method, df = NA)
                }
                lavaSearch2::getVarCov2(model, data = data, individual = x, plot = FALSE)
            })
            V <- Matrix::bdiag(ls.V)             
        }
        iV <- solve(V)
    }
    
    ## R2
    Y <- data[[name.endo]]
    value0 <- mean(Y)    
    value <- predict(model, type = "response")
    R2 <- .calcR2(residuals = Y - value, residuals0 = Y - value0, iV = iV)

    ## McFadden R2
    iff <- stats::update(ff, paste0(".~1"))
    imodel <- try(stats::update(model, iff, data = data), silent = TRUE)
    if(class(imodel)!="try-error"){
        LR <- 2 * as.numeric(stats::logLik(model)-stats::logLik(imodel))
        R2.LR <- 1-exp(-LR/n)
    }else{
        R2.LR <- NA
    }
    
    ## partial R2
    pR2 <- setNames(numeric(n.exo),name.exo)
    pR2.LR <- setNames(numeric(n.exo),name.exo)
    for(iX in 1:n.exo){ # iX <- 1
        if(trace){ cat(iX,"(",name.exo[iX],") ", sep = "") }
        iff <- stats::update(ff, paste0(".~.-",name.exo[iX]))
        imodel <- try(stats::update(model, iff, data = data), silent = TRUE)
        if(class(imodel)!="try-error"){
            value0 <- predict(imodel, type = "response")
            pR2[iX] <- .calcR2(residuals = Y - value, residuals0 = Y - value0, iV = iV)
            LR <- 2 * as.numeric(stats::logLik(model)-stats::logLik(imodel))
            pR2.LR[iX] <- 1-exp(-LR/n)
        }
    }
    if(trace){ cat("\n") }

    return(list(R2.Buse = as.numeric(R2),
                R2.McFadden = R2.LR,
                partialR2.Buse = pR2,
                partialR2.McFadden = pR2.LR))
}

## * FCT.calcR2
.calcR2 <- function(residuals, residuals0, iV){
    if(is.null(iV)){
        SS.residual <- sum(residuals^2)
        SS.total <- sum(residuals0^2)
    }else{
        SS.residual <- t(residuals) %*% iV %*% residuals
        SS.total <- t(residuals0) %*% iV %*% residuals0
    }
    R2 <- 1 - SS.residual/SS.total
    return(R2)
}

#----------------------------------------------------------------------
### R2.R ends here




