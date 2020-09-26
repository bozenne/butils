### glmEff.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 26 2020 (18:00) 
## Version: 
## Last-Updated: sep 26 2020 (19:22) 
##           By: Brice Ozenne
##     Update #: 11
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
#'
#' library(lava)
#'
#' 
#' @export
glmEff <- function(outcome, treatment, covariates, data, family.outcome){
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

    family.outcome <- match.arg(family.outcome, c("gaussian","binomial"))
    
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
    
    ls.nuisance <- lapply(1:length(level.treatment), function(iT){
        glm(ff.nuisance, data = data[ls.index.treatment[[iT]],], family = family.outcome)
    })
    ls.fitted.nuisance <- lapply(ls.nuisance, predict, type = "response")
    M.fitted.nuisance <- cbind(1-unlist(ls.fitted.nuisance)[order(unlist(ls.index.treatment))],
                               unlist(ls.fitted.nuisance)[order(unlist(ls.index.treatment))])
    
    ## ** augmented estimating equation
    pi <- sapply(level.treatment, function(iT){ mean(data[[treatment]]==iT)})
    M.Impi <- do.call(cbind,lapply(1:length(level.treatment), function(iT){ (data[[treatment]]==level.treatment[iT]) - pi[iT]  }))
    ls.X <- lapply(level.treatment, function(iT){
        iData <- copy(data)
        iData[[treatment]] <- iT
        X <- model.matrix(ff.treatment, iData)
    })

    browser()
  
    if(family.outcome == "gaussian"){
        e0 <- glm(ff.treatment, data = data, family = family.outcome)

        a <- mean(data$Y[data$T==0]) + (1/sum(data$T==0)) * sum((data$T-mean(data$T==1))*M.fitted.nuisance[,1])
        b <- mean(data$Y[data$T==1]) - (1/sum(data$T==1)) * sum((data$T-mean(data$T==1))*M.fitted.nuisance[,2])
        b-a
        mean(data$Y[data$T==0]) - (1/sum(data$T==0)) * sum(M.Impi[,1]*M.fitted.nuisance[,1])
        mean(data$Y[data$T==1]) - (1/sum(data$T==1)) * sum(M.Impi[,2]*M.fitted.nuisance[,2])

        coef(e0) - colSums(M.Impi*M.fitted.nuisance)/(pi*n.obs)
        ## rhs <- t(X) %*% Y - Reduce("+",lapply(1:length(level.treatment), function(iT){ t(ls.X[[iT]]) %*% (M.fitted.nuisance[,iT] * M.Impi[,iT]) }))
        ## lhs <- t(X) %*% X - Reduce("+",lapply(1:length(level.treatment), function(iT){ t(ls.X[[iT]]) %*% sweep(ls.X[[iT]], MARGIN = 1, FUN = "*", STATS = M.Impi[,iT])}))
        ## out <- (solve(lhs) %*% rhs)[,1]

        
        ## alternative computation
        
        d1 <- function(theta){ ## theta <- coef(e0)
            term.score <- - t(X) %*% (Y - X %*% theta)
            term.aug <- - t(ls.X[[2]]) %*% ((M.fitted.nuisance[,2] - ls.X[[2]] %*% theta) * M.Impi[,2]) + t(ls.X[[1]]) %*% ((M.fitted.nuisance[,1] - ls.X[[1]] %*% theta) * M.Impi[,2])
            ## term.aug <- Reduce("+",lapply(1:length(level.treatment), function(iT){ ## iT <- 2
            ##     - t(ls.X[[iT]]) %*% ((M.fitted.nuisance[,iT] - ls.X[[iT]] %*% theta) * M.Impi[,iT])
            ## }))
            return(term.score - term.aug)
        }
        ## d2 <- function(theta){             ## theta <- coef(e0)
        ##     term.score <- t(X) %*% X
        ##     term.aug <- Reduce("+",lapply(1:length(level.treatment), function(iT){ ## iT <- 1
        ##         t(ls.X[[iT]]) %*% sweep(ls.X[[iT]], MARGIN = 1, FUN = "*", STATS = M.Impi[,iT])
        ##     }))
        ##     return(term.score + term.aug)
        ## }

        ## e <- rootSolve::multiroot(f = d1, jacfunc = d2, start = coef(e0))
        e <- rootSolve::multiroot(f = d1, start = coef(e0))
        out <- e$root

        ## t(X) %*% (Y - X %*% theta) - t(ls.X[[1]]) %*% ((M.fitted.nuisance[,1] - ls.X[[1]] %*% theta) * M.Impi[,1]) - t(ls.X[[2]]) %*% ((M.fitted.nuisance[,2] - ls.X[[2]] %*% theta) * M.Impi[,2])
        

    }else if(family.outcome == "binomial"){
        ## t(X) %*% (Y - 1/(1+exp(-X %*% theta)))
        e0 <- glm(ff.treatment, data = data, family = family.outcome)
        
        d1 <- function(theta){
            term.score <- - t(X) %*% (Y - 1/(1+exp(-X %*% theta)))
            term.aug <- Reduce("+",lapply(1:length(level.treatment), function(iT){ ## iT <- 2
                - t(ls.X[[iT]]) %*% ((M.fitted.nuisance[,iT] - 1/(1+exp(-ls.X[[iT]] %*% theta))) * M.Impi[,iT])
            }))
            return(term.score + term.aug)
        }
        
        d2 <- function(theta){             ## theta <- coef(e0)
            term.score <- t(X) %*% sweep(X, MARGIN = 1, FUN = "*", STATS = exp(-X %*% theta)/(1+exp(-X %*% theta))^2)
            term.aug <- Reduce("+",lapply(1:length(level.treatment), function(iT){ ## iT <- 1
                ls.X[[iT]] %*% sweep(ls.X[[iT]], MARGIN = 1, FUN = "*",
                                     STATS = exp(-ls.X[[iT]] %*% theta) * M.Impi[,iT]/(1+exp(-ls.X[[iT]] %*% theta))^2)
            }))
            return(term.score + term.aug)
        }

        e <- rootSolve::multiroot(f = d1, jacfunc = d2, start = coef(e0))
        ## e <- rootSolve::multiroot(f = d1, start = coef(e0))
        out <- e$root
    }
    return(out)
}


######################################################################
### glmEff.R ends here
