### glmEff.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 26 2020 (18:00) 
## Version: 
## Last-Updated: sep 28 2020 (18:23) 
##           By: Brice Ozenne
##     Update #: 37
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * glmEff (documentation)
#' @title Efficient Group Comparison Using Baseline Adjustment
#' @description Efficient group comparison in a randomized experiment using baseline adjustment.
#' @name glmEff
#' 
#' @param outcome [character] the name of the outcome variable in the dataset
#' @param treatment [character] the name of the treatment variable in the dataset
#' @param covariates [character] the name of the baseline variables in the dataset
#' @param data [data.frame] the dataset
#' @param family.outcome [character] the type of outcome: can be continuous (\code{"gaussian"}) or binary (\code{"binomial"}).
#'
#' @references Zhang et al. 2010. Increasing the Efficiency of Prevention Trials by Incorporating Baseline Covariates. Statistical Communications in Infectious Diseases

## * glmEff (examples)
#' @rdname glmEff
#' @examples
#' library(lava)
#'
#' ## continuous case
#' m.cont <- lvm(Y~T+3*X1-3*X2)
#' distribution(m.cont,~T) <- binomial.lvm(size=1,p=0.5)
#' distribution(m.cont,~X1) <- binomial.lvm(size=1,p=0.25)
#'
#' set.seed(10)
#' d <- lava::sim(m.cont, n = 100)
#'
#' m.glmeff <- glmEff(outcome = "Y", treatment = "T", covariates = c("X1","X2"), data = d,
#'                    family.outcome = "gaussian")
#' m.glmeff ## efficient
#' tapply(d$Y,d$T,mean) ## naive
#'
#' ## binary case
#' m.bin <- lvm(Y~T+3*X1-3*X2)
#' distribution(m.bin,~Y) <- binomial.lvm(size=1,p=0.5)
#' distribution(m.bin,~T) <- binomial.lvm(size=1,p=0.5)
#' distribution(m.bin,~X1) <- binomial.lvm(size=1,p=0.25)
#'
#' set.seed(10)
#' d <- lava::sim(m.bin, n = 1000)
#'
#' m.glmeff <- glmEff(outcome = "Y", treatment = "T", covariates = c("X1","X2"), data = d,
#'                    family.outcome = "binomial")
#' m.glmeff ## efficient
#' tapply(d$Y,d$T,mean) ## naive

## * glmEff (code)
#' @export
glmEff <- function(outcome, treatment, covariates, data, family.outcome,
                   method = "explicit"){
    require(data.table)
    
    ## ** check arguments
    data <- as.data.frame(data)
    if(length(outcome)!=1){
        stop("Argument \'outcome\' must have length 1 \n")
    }
    if(outcome %in% names(data) == FALSE){
        stop("Argument \'outcome\' not found in argument \'data\' \n")
    }
    if(length(treatment)!=1){
        stop("Argument \'treatment\' must have length 1 \n")
    }
    if(treatment %in% names(data) == FALSE){
        stop("Argument \'treatment\' not found in argument \'data\' \n")
    }
    if(length(covariates)==0){
        stop("Argument \'covariates\' must have length at least 1 \n")
    }
    if(any(covariates %in% names(data) == FALSE)){
        stop("Argument \'covariates\' not found in argument \'data\' \n")
    }
    if(any(is.na(data))){
        stop("Do not currently handle missing data \n")
    }
    if("offset.aug" %in% names(data)){
        stop("Argument \'data\' should not contain a column called \"offset.aug\". \n",
             "This name is used internally.")
    }
    if("intercept" %in% names(data)){
        stop("Argument \'data\' should not contain a column called \"intercept\". \n",
             "This name is used internally.")
    }
    method <- match.arg(method, c("explicit","multiroot"))

    family.outcome <- match.arg(family.outcome, c("gaussian","binomial"))
    link.outcome <- switch(family.outcome,
                           "gaussian" = function(x){x},
                           "binomial" = function(x){1/(1+exp(-x))})
    ilink.outcome <- switch(family.outcome,
                            "gaussian" = function(x){x},
                            "binomial" = function(x){log(x/(1-x))})
    
    ## ** prepare
    level.treatment <- sort(unique(data[[treatment]]))
    ls.index.treatment <- setNames(lapply(level.treatment, function(iT){
        which(data[[treatment]]==iT)
    }),level.treatment)
    n.obs <- NROW(data)
    
    ## ** estimating equation
    ff.treatment <- as.formula(paste0(outcome,"~",treatment))
    X <- model.matrix(ff.treatment, data)
    Y <- data[[outcome]]
    
    ## ** fit nuisance model
    ff.nuisance <- as.formula(paste0(outcome,"~",paste0(covariates,collapse="+")))
    
    M.nuisance <- do.call(cbind,lapply(1:length(level.treatment), function(iT){
        iGLM <- glm(ff.nuisance, data = data[ls.index.treatment[[iT]],], family = family.outcome)
        predict(iGLM, newdata = data, type = "response")
    }))
    
    ## ** augmented estimating equation
    pi <- sapply(level.treatment, function(iT){ mean(data[[treatment]]==iT)})
    M.Impi <- do.call(cbind,lapply(1:length(level.treatment), function(iT){ (data[[treatment]]==level.treatment[iT]) - pi[iT]  }))
    ls.X <- lapply(level.treatment, function(iT){
        iData <- copy(data)
        iData[[treatment]] <- iT
        X <- model.matrix(ff.treatment, iData)
    })

    e0 <- glm(ff.treatment, data = data, family = family.outcome)
    
    if(method == "multiroot"){

        d1 <- function(theta){ ## theta <- coef(e0)
            term.score <- - t(X) %*% (Y - link.outcome(X %*% theta))
            eXb <- link.outcome(do.call(cbind,lapply(ls.X,`%*%`,theta)))
            h <- (M.nuisance - eXb) * M.Impi
            term.aug <- Reduce("+",lapply(1:length(level.treatment), function(iT){- t(ls.X[[iT]]) %*% h[,iT]}))
            return(term.score - term.aug)
        }
        e <- rootSolve::multiroot(f = d1, start = coef(e0))
        out <- setNames(c(cumsum(e$root),utils::tail(e$root,1)),c(level.treatment, "difference"))

        ## Gaussian case
        ## d2 <- function(theta){             ## theta <- coef(e0)
        ##     term.score <- t(X) %*% X
        ##     term.aug <- Reduce("+",lapply(1:length(level.treatment), function(iT){ ## iT <- 1
        ##         t(ls.X[[iT]]) %*% sweep(ls.X[[iT]], MARGIN = 1, FUN = "*", STATS = M.Impi[,iT])
        ##     }))
        ##     return(term.score + term.aug)
        ## }

        ## e <- rootSolve::multiroot(f = d1, jacfunc = d2, start = coef(e0))

    }else{
        e <-  link.outcome(cumsum(coef(e0))) - 1/(pi*n.obs) * colSums(M.Impi*M.nuisance)
        out <- setNames(c(e, diff(e)), c(level.treatment, "difference"))
    }        
    return(out)
}


######################################################################
### glmEff.R ends here
