### pairedCompMean.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep  6 2021 (16:56) 
## Version: 
## Last-Updated: sep  7 2021 (11:33) 
##           By: Brice Ozenne
##     Update #: 91
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
##' @param col.factor [character] column in the dataset indicating the factor variable, to stratify the comparison on.
##' @param col.value [character] column in the dataset indicating the value of the measurement.
##' @param col.id [character] column in the dataset indicating the patient identity.
##' @param n.perm [numeric] number of permutation to be performed.
##' @param seed [integer,>0] set initial state of the random number generation (passed to \code{set.seed}).
##' @param trace [logical] should a progress bar be display displaying the execution of the permutations.
##'
##' @details Single step max test adjustment is performed to adjust for multiple comparisons. It account for correlation between test statistics.
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
##'                col.treatment = "treatment", col.factor = "score", col.value = "value", col.cluster = "id",
##'                seed = NULL)
##' resPerm
##'



## * pairedCompMean (code)
##' @rdname pairedCompMean
##' @export
pairedCompMean <- function(data, col.treatment, col.factor, col.value, col.cluster,
                           n.perm = 1e3, method.adj = "step-down", seed = NULL, trace = TRUE){

    if(!is.null(seed)){set.seed(seed)}
    Utreatment <- unique(data[[col.treatment]])
    n.treatment <- length(Utreatment)
    Ufactor <- unique(data[[col.factor]])
    data <- as.data.frame(data)
    data <- data[order(data[[col.cluster]],data[[col.factor]]),,drop=FALSE]
    
    grid <- expand.grid(Utreatment[-1],
                        Ufactor,
                        stringsAsFactors = FALSE)
    colnames(grid) <- c(col.treatment,col.factor)
    n.grid <- NROW(grid)

    method.adj <- match.arg(method.adj, c("step-down","single-step"))
    
    ## ** point estimate
    ls.ttest <- lapply(1:n.grid, function(iG){ ## iG <- 1
        .pairedttest(x = data[data[[col.treatment]]==grid[iG,col.treatment] & data[[col.factor]] == grid[iG,col.factor],col.value],
                     y = data[data[[col.treatment]]==Utreatment[1] & data[[col.factor]] == grid[iG,col.factor],col.value])
    })
    df.ttest <- data.frame(grid, ref = Utreatment[1], do.call(rbind,ls.ttest))

    ## ** permutation
    ls.indexID <- tapply(1:NROW(data),data[[col.cluster]],function(iVec){iVec})
    n.id <- length(ls.indexID)
    ls.indexFactor <- tapply(1:NROW(data),data[[col.factor]],function(iVec){iVec})

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
        
        ## ttest
        iLS.ttest <- lapply(1:n.grid, function(iG){ ## iG <- 1
            .pairedttest(value[intersect(iLS.indexTreatment[[grid[iG,col.treatment]]],ls.indexFactor[[grid[iG,col.factor]]])],
                         value[intersect(iLS.indexTreatment[[Utreatment[1]]],ls.indexFactor[[grid[iG,col.factor]]])])
        })
        iDF.ttest <- data.frame(grid, ref = Utreatment[1], do.call(rbind,iLS.ttest))
        if(trace){
            utils::setTxtProgressBar(pb, iPerm)
        }
        ls.permttest[[iPerm]] <- cbind(perm = iPerm,iDF.ttest)
    }
    if(trace){
        close(pb)
    }
    
    ## ** post-process
    M.permttest <- do.call(rbind, ls.permttest)

    df.ttest$perm.p.value <- NA
    for(iGrid in 1:n.grid){ ## iGrid <- 1
        vec.H0 <- abs(M.permttest[M.permttest[[col.treatment]] == grid[iGrid,col.treatment] & M.permttest[[col.factor]] == grid[iGrid,col.factor],"statistic"])
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
        
        for(iI in 1:length(iIndex.treatment2)){ ## iI <- 3
            iII <- iIndex.treatment2[iI]
            if(method.adj=="step-down"){
                iM.permttest.remain <- iM.permttest[iM.permttest$score %in% df.ttest[iIndex.treatment2[iI:length(iIndex.treatment2)],col.factor],,drop=FALSE]
                vec.H0 <- tapply(abs(iM.permttest.remain$statistic),iM.permttest.remain$perm, max)
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
