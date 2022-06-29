### pairedCompMean.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep  6 2021 (16:56) 
## Version: 
## Last-Updated: jun 29 2022 (15:24) 
##           By: Brice Ozenne
##     Update #: 110
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * pairedCompMean (documentation)
##' @title Permutation Test for the Mean with Paired Data
##' @description Permutation test for comparing mean between within cluster with efficient adjustment for multiple comparisons
##' @name pairedCompMean
##' @param data [data.frame] dataset with in row the samples.
##' @param col.treatment [character] column in the dataset indicating the treatment variable, that should be permuted within-individuals.
##' @param col.strata [character] column in the dataset indicating the factor variable, to stratify the comparison on.
##' @param col.value [character] column in the dataset indicating the value of the measurement.
##' @param col.cluster [character] column in the dataset indicating the cluster level (e.g. patient identity).
##' @param n.perm [numeric] number of permutation to be performed. 
##' @param method.adj [character] either \code{"single-step"} or \code{"step-down"} max adjustment.
##' Both order the hypothesis to test based on the evidence shown by the data but the latter provides extra statistical power for rejecting the 2nd to last hypotheses.
##' @param seed [integer,>0] set initial state of the random number generation (passed to \code{set.seed}).
##' @param trace [logical] should a progress bar be display displaying the execution of the permutations.
##' 
##' @examples
##' library(mvtnorm)
##' library(reshape2)
##'
##' ## generate data
##' n <- 25
##' 
##' Sigma1 <- (diag(1-0.5,5,5)+0.5) * tcrossprod(c(1,2,1,2,3))
##' Sigma2 <- (diag(1-0.5,5,5)+0.5) * tcrossprod(c(1,2,1,2,3))
##' Sigma3 <- (diag(1-0.5,5,5)+0.5) * tcrossprod(c(1,2,1,2,3))
##' Sigma <- as.matrix(Matrix::bdiag(Sigma1,Sigma2,Sigma3))+0.25
##' mu1 <- c(1,2,1,2,1.5)
##' mu <- c(mu1, mu1+1.5,mu1-1)
##'
##' set.seed(10)
##' dfW <- data.frame(1:n,mvtnorm::rmvnorm(n, mean = mu, sigma = Sigma))
##' colnames(dfW) <- c("id",paste0("score",1:5,"_B"),
##'                    paste0("score",1:5,"_P"),
##'                    paste0("score",1:5,"_K"))
##' dfL <- melt(dfW, id.vars = "id")
##' dfL$score <- sapply(strsplit(x=as.character(dfL$variable),split="_",fixed = TRUE),"[",1)
##' dfL$treatment <- sapply(strsplit(x=as.character(dfL$variable),split="_",fixed = TRUE),"[",2)
##'
##' ## test mean
##' t.test(dfW$score1_P,dfW$score1_B,paired = TRUE)
##' t.test(dfW$score5_P,dfW$score5_B,paired = TRUE)
##'
##' if(require(MKinfer)){
##' perm.t.test(dfW$score1_P,dfW$score1_B,paired = TRUE)
##' perm.t.test(dfW$score2_P,dfW$score2_B,paired = TRUE)
##' perm.t.test(dfW$score3_P,dfW$score3_B,paired = TRUE)
##' perm.t.test(dfW$score4_P,dfW$score4_B,paired = TRUE)
##' perm.t.test(dfW$score5_P,dfW$score5_B,paired = TRUE)
##' }
##' 
##' resPerm <- pairedCompMean(dfL[dfL$treatment %in% c("P","B"),], n.perm = 1e3,
##'                col.treatment = "treatment", col.strata = "score", col.value = "value", col.cluster = "id",
##'                seed = NULL)
##' resPerm
##'



## * pairedCompMean (code)
##' @rdname pairedCompMean
##' @export
pairedCompMean <- function(data, col.treatment, col.strata, col.value, col.cluster,
                           n.perm = 1e4, method.adj = "single-step", seed = NULL, trace = TRUE){


    ## ** normalize user input
    if(!is.null(seed)){set.seed(seed)}
    if(!inherits(data,"data.frame")){
        stop("Argument \'data\' must be or inherit from data.frame.")
    }
    if(col.treatment %in% names(data) == FALSE){
        stop("Argument \'col.treatment\' must correspond to a column in argument \'data\'.")
    }
    if(col.strata %in% names(data) == FALSE){
        stop("Argument \'col.strata\' must correspond to a column in argument \'data\'.")
    }
    if(col.value %in% names(data) == FALSE){
        stop("Argument \'col.value\' must correspond to a column in argument \'data\'.")
    }
    if(col.cluster %in% names(data) == FALSE){
        stop("Argument \'col.cluster\' must correspond to a column in argument \'data\'.")
    }

    ## ** prepare
    Utreatment <- unique(data[[col.treatment]])
    n.treatment <- length(Utreatment)
    Ufactor <- unique(data[[col.strata]])
    data <- as.data.frame(data)
    data <- data[order(data[[col.cluster]],data[[col.strata]]),,drop=FALSE]

    grid <- expand.grid(Utreatment[-1],
                        Ufactor,
                        stringsAsFactors = FALSE)
    colnames(grid) <- c(col.treatment,col.strata)
    n.grid <- NROW(grid)

    method.adj <- match.arg(method.adj, c("step-down","single-step"))
    
    ## ** point estimate
    ls.ttest <- lapply(1:n.grid, function(iG){ ## iG <- 1
        .pairedttest(x = data[data[[col.treatment]]==grid[iG,col.treatment] & data[[col.strata]] == grid[iG,col.strata],col.value],
                     y = data[data[[col.treatment]]==Utreatment[1] & data[[col.strata]] == grid[iG,col.strata],col.value])
    })
    df.ttest <- data.frame(grid, ref = Utreatment[1], do.call(rbind,ls.ttest))

    ## ** permutation
    ls.indexID <- tapply(1:NROW(data),data[[col.cluster]],function(iVec){iVec})
    n.id <- length(ls.indexID)
    ls.indexFactor <- tapply(1:NROW(data),data[[col.strata]],function(iVec){iVec})

    allPerm <- .allPermutation(Utreatment)
    colnames(allPerm) <- Utreatment
    n.allPerm <- NROW(allPerm)
    
    ls.permttest <- vector(mode = "list", length = n.perm)
    if(trace){
        pb <- utils::txtProgressBar(max = n.perm)
    }
    value <- data[[col.value]]
    
    for(iPerm in 1:n.perm){
        iTreatment <- data[[col.treatment]]
        iIndexPerm <- sample.int(n.allPerm, size = n.id, replace = TRUE)

        ## permute treatment within individual
        for(iId in 1:n.id){ ## iId <- 1
            iTreatment[ls.indexID[[iId]]] <- allPerm[iIndexPerm[iId],][iTreatment[ls.indexID[[iId]]]]
            ## microbenchmark::microbenchmark(a = unname(setNames(sample(Utreatment),Utreatment)[iData[ls.indexID[[iId]],col.treatment]]),
            ##                                b = as.character(factor(iData[ls.indexID[[iId]],col.treatment], levels = Utreatment, labels = sample(Utreatment))))
        }
        iLS.indexTreatment <- tapply(1:NROW(data),iTreatment,function(iVec){iVec})
        iLS.valueRef <- stats::setNames(lapply(Ufactor, function(iFactor){
            value[intersect(iLS.indexTreatment[[Utreatment[1]]],ls.indexFactor[[iFactor]])]
        }), Ufactor)
        

        ## ttest
        iLS.ttest <- lapply(1:n.grid, function(iG){ ## iG <- 1
            .pairedttest(value[intersect(iLS.indexTreatment[[grid[iG,col.treatment]]],ls.indexFactor[[grid[iG,col.strata]]])],
                         iLS.valueRef[[grid[iG,col.strata]]])
        })
        ls.permttest[[iPerm]] <- data.frame(grid, perm = iPerm, ref = Utreatment[1], do.call(rbind,iLS.ttest))
        if(trace){
            utils::setTxtProgressBar(pb, iPerm)
        }
    }
    if(trace){
        close(pb)
    }
    
    ## ** post-process
    M.permttest <- do.call(rbind, ls.permttest)

    df.ttest$perm.p.value <- NA
    for(iGrid in 1:n.grid){ ## iGrid <- 1
        vec.H0 <- abs(M.permttest[M.permttest[[col.treatment]] == grid[iGrid,col.treatment] & M.permttest[[col.strata]] == grid[iGrid,col.strata],"statistic"])
        vec.H1 <- abs(df.ttest[iGrid,"statistic"])
        df.ttest[iGrid,"perm.p.value"] <- (sum(vec.H0>vec.H1)+1)/(n.perm+1)
    }
    df.ttest$adj.p.value <- NA
    for(iTreatment in 2:n.treatment){ ## iTreatment <- 2

        iIndex.treatment <- which(df.ttest[[col.treatment]]==Utreatment[iTreatment])
        iIndex.treatment2 <- iIndex.treatment[order(abs(df.ttest[iIndex.treatment,"statistic"]), decreasing = TRUE)]
        
        iM.permttest <- M.permttest[M.permttest[[col.treatment]] == Utreatment[iTreatment],]

        if(method.adj=="single-step"){
            vec.H0 <- tapply(abs(iM.permttest$statistic),iM.permttest$perm, max)
        }
        
        for(iI in 1:length(iIndex.treatment2)){ ## iI <- 1
            iII <- iIndex.treatment2[iI]
            if(method.adj=="step-down"){
                iM.permttest.remain <- iM.permttest[iM.permttest[[col.strata]] %in% df.ttest[iIndex.treatment2[iI:length(iIndex.treatment2)],col.strata],,drop=FALSE]
                vec.H0 <- tapply(abs(iM.permttest.remain$statistic),iM.permttest.remain$perm, max)
                ## vec.H0 - tapply(abs(iM.permttest.remain$statistic),iM.permttest.remain$perm, max)
            }
            df.ttest[iII,"adj.p.value"] <- max((sum(vec.H0>abs(df.ttest[iII,"statistic"]))+1)/(n.perm+1),
                                               df.ttest[iIndex.treatment2,"adj.p.value"], na.rm = TRUE)
        }
        
    }    

    ## ** export
    return(df.ttest)
}

## * .pairedttest
##' @title Fast Paired t-test
##' @description Fast paired t-test
##' @param x vector of numeric values 
##' @param y vector of numeric values (reference)
##' @examples
##' X <- rnorm(1e2)
##' Y <- rnorm(1e2)
##' t.test(X,Y,paired = TRUE)
##' .pairedttest(X,Y)
##' microbenchmark::microbenchmark(t.test(X,Y,paired = TRUE), .pairedttest(X,Y))
##' @keywords internal
.pairedttest <- function(x,y){
    diff <- na.omit(x - y)
    n.obs <- length(diff)
    meanDiff <- sum(diff)/n.obs
    sdDiff <- sqrt(sum((diff-meanDiff)^2)/(n.obs-1))
    statistic <- sqrt(n.obs)*meanDiff/sdDiff
    pvalue <- 2 * stats::pt(abs(statistic), n.obs - 1, lower.tail = FALSE)
    out <- c(mean.x = mean(x,na.rm = TRUE),
             mean.y = mean(y,na.rm = TRUE),
             estimate =  meanDiff,
             statistic = statistic,
             asym.p.value = pvalue)
    return(out)
}

## * .allPermutation
##' @title Generate all possible permutations
##' @description Generate all possible permutations
##' @param x vector of values
##' @keywords internal
##' @references https://stackoverflow.com/questions/11095992/generating-all-distinct-permutations-of-a-list-in-r
##' @examples
##' .allPermutation(1)
##' .allPermutation(1:2)
##' .allPermutation(1:3)
.allPermutation <- function(x) {
    if (length(x) == 1) {
        return(x)
    }
    else {
        res <- matrix(nrow = 0, ncol = length(x))
        for (i in seq_along(x)) {
            res <- rbind(res, cbind(x[i], Recall(x[-i])))
        }
        return(res)
    }
}
##----------------------------------------------------------------------
### pairedCompMean.R ends here
