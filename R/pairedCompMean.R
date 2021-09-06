### pairedCompMean.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep  6 2021 (16:56) 
## Version: 
## Last-Updated: sep  6 2021 (17:12) 
##           By: Brice Ozenne
##     Update #: 10
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
##' @param col.group [character] column in the dataset indicating the treatment variable.
##' @param col.factor [character] column in the dataset indicating the factor variable, to stratify the comparison on.
##' @param col.value [character] column in the dataset indicating the value of the measurement.
##' @param col.id [character] column in the dataset indicating the patient identity.
##' @param col.id [character] column in the dataset indicating the patient identity.
##' 
##' @examples
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
##' dfW <- data.frame(1:n,rmvnorm(n, mean = mu, sigma = Sigma))
##' colnames(dfW) <- c("id",paste0("score",1:5,"_B"),
##'                    paste0("score",1:5,"_P"),
##'                    paste0("score",1:5,"_K"))
##' dfL <- melt(dfW, id.vars = "id")
##' dfL$score <- sapply(strsplit(x=as.character(dfL$variable),split="_",fixed = TRUE),"[",1)
##' dfL$treatment <- sapply(strsplit(x=as.character(dfL$variable),split="_",fixed = TRUE),"[",2)
##'
##' ## test mean
##' pairedCompMean(dfL[dfL$treatment %in% c("P","B"),])
##'



## * pairedCompMean (code)
##' @rdname pairedCompMean
##' @export
pairedCompMean <- function(data, col.group = "treatment", col.factor = "score", col.value = "value", col.cluster = "id",
                     n.perm = 1e3, trace = TRUE){
    Ugroup <- unique(data[[col.group]])
    n.group <- length(Ugroup)
    Ufactor <- unique(data[[col.factor]])
    data <- as.data.frame(data)
    data <- data[order(data[[col.cluster]],data[[col.factor]]),,drop=FALSE]
    
    grid <- expand.grid(group = Ugroup[-1],
                        factor = Ufactor)
    n.grid <- NROW(grid)

    ## ** point estimate
    ls.ttest <- lapply(1:n.grid, function(iG){ ## iG <- 1
        iTT <- t.test(x = data[data[[col.group]]==grid$group[iG] & data[[col.factor]] == grid$factor[iG],col.value],
                      y = data[data[[col.group]]==Ugroup[1] & data[[col.factor]] == grid$factor[iG],col.value])
        c(mean_group = as.double(iTT$estimate[1]),
          mean_ref = as.double(iTT$estimate[2]),
          estimate = as.double(-diff(iTT$estimate)),
          statistic = as.double(iTT$statistic),
          lower = iTT$conf.int[1],
          upper = iTT$conf.int[2],
          p.value = iTT$p.value)
    })
    df.ttest <- data.frame(grid, ref = Ugroup[1], do.call(rbind,ls.ttest))

    ## ** permutation
    iOutPerm <- NULL
    if(trace){
        pb <- utils::txtProgressBar(max = n.perm)
    }
    for(iPerm in 1:n.perm){
        iData <- data
        ## permute group within individual
        iData[[col.group]] <- unlist(
            tapply(iData[[col.group]],iData[[col.cluster]],function(iG){
                as.character(factor(iG,levels = Ugroup, labels = sample(Ugroup)))
            })
        )
        ## ttest
        iLS.ttest <- lapply(1:n.grid, function(iG){ ## iG <- 1
            iTT <- t.test(x = iData[iData[[col.group]]==grid$group[iG] & data[[col.factor]] == grid$factor[iG],col.value],
                          y = iData[iData[[col.group]]==Ugroup[1] & data[[col.factor]] == grid$factor[iG],col.value])
            c(mean_group = as.double(iTT$estimate[1]),
              mean_ref = as.double(iTT$estimate[2]),
              estimate = as.double(-diff(iTT$estimate)),
              statistic = as.double(iTT$statistic),
              lower = iTT$conf.int[1],
              upper = iTT$conf.int[2],
              asym.p.value = iTT$p.value)
        })
        iDF.ttest <- data.frame(grid, ref = Ugroup[1], do.call(rbind,iLS.ttest))
        if(trace){
            utils::setTxtProgressBar(pb, iPerm)
        }
        iOutPerm <- rbind(iOutPerm,
                          cbind(perm = iPerm,iDF.ttest))
    }
    if(trace){
        close(pb)
    }

    ## ** post-process
    df.ttest$perm.p.value <- NA
    for(iGrid in 1:n.grid){ ## iGrid <- 1
        vec.H0 <- abs(iOutPerm[iOutPerm[["factor"]] == grid[iGrid,"factor"] & iOutPerm[["group"]] == grid[iGrid,"group"],"statistic"])
        vec.H1 <- abs(df.ttest[iGrid,"statistic"])
        df.ttest[iGrid,"perm.p.value"] <- mean(vec.H0>vec.H1)
    }
    df.ttest$adj.p.value <- NA
    for(iGroup in 2:n.group){ ## iGroup <- 2
        iiOutPerm <- iOutPerm[iOutPerm$group == Ugroup[iGroup],]
        vec.H0 <- tapply(abs(iiOutPerm$statistic),iiOutPerm$perm, max)
        df.ttest[df.ttest$group==Ugroup[iGroup],"adj.p.value"] <- sapply(abs(df.ttest[df.ttest$group==Ugroup[iGroup],"statistic"]), function(iT){mean(vec.H0>iT)})
    }    

    ## ** export
    return(df.ttest)
}


##----------------------------------------------------------------------
### pairedCompMean.R ends here
