### compMean.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 10 2020 (10:27) 
## Version: 
## Last-Updated: okt  3 2023 (10:32) 
##           By: Brice Ozenne
##     Update #: 62
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * compMean (documentation)
##' @title Permutation Test for the Mean
##' @description Permutation test for comparing mean between groups with efficient adjustment for multiple comparisons
##' @param Y [matrix] Matrix where each column correspond to a difference outcome.
##' @param group [character vector] Group variable.
##' @param time [character vector] Optional time variable.
##' Can only take two different values. The mean difference and the mean change over time between the two groups will then be tested.
##' @param permutation.type [character] Should the group labels (\code{"group"}) or the outcome labels (\code{"Y"}) be permuted.
##' In the latter case, exactly two outcomes should be specified and the function will test whether the group differences are the same for both outcomes.
##' @param n.perm [integer, >0] number of permutations.
##' @param na.rm [logical] Should missing values be ignored when evaluating mean and variance (complete case analysis).
##' @param cl [integer, >0] Optional cluster or number of cores used to run the permutations in parallel.
##' No available for Windows OS.
##' @param trace [logical] Should a progress bar be displayed to follow the permutations.
##'
##' @details Several adjustment for multiple comparisons are available via the print function \itemize{
##' \item Holm: only accounts for the ordering of the tests to gain efficiency compared to a Bonferroni adjustment.
##' \item max: only accounts for the correlation between the tests to gain efficiency compared to a Bonferroni adjustment.
##' \item Step-down max: accounts both for the correlation between the tests and their ordering to gain efficiency compared to a Bonferroni adjustment.
##' }
##' 
##' @references
##' Westfall, P. H. and  Troendle, J. F. Multiple Testing with Minimal Assumptions, Biometrical Journal 50 (2008) 5, 745–755 DOI: 10.1002/bimj.200710456b.
##' Dudoit S., Shaffer J. P. and Boldrick J. C. Multiple  Hypothesis Testing  in Microarray  Experiments. Statistical Science 2003, Vol. 18, No. 1, 71–103.

## * compMean (examples)
#' @examples
#' library(data.table)
#' set.seed(11)
#' n <- 100
#'
#' #### 1- H0 under the global null ####
#' dtH0 <- rbind(data.table(Y1 = rnorm(n), Y2 = rnorm(n), group = "C", time = "D1"),
#'               data.table(Y1 = rnorm(n), Y2 = rnorm(n), group = "T", time = "D1"),
#'               data.table(Y1 = rnorm(n), Y2 = rnorm(n), group = "C", time = "W1"),
#'               data.table(Y1 = rnorm(n), Y2 = rnorm(n), group = "T", time = "W1"))
#'
#' dtH0.mean <- dtH0[,.(meanY1 = mean(Y1)),by=c("group","time")]
#' dtH0.mean[group=="T", diff(meanY1)] - dtH0.mean[group=="C", diff(meanY1)]
#' 
#' res1 <- compMean(Y = cbind(Y1 = dtH0$Y1, Y2 = dtH0$Y2),
#'                  group = dtH0$group, time = dtH0$time, permutation.type = "group")
#' print(res1)
#' print(res1, method = "max")
#' print(res1, method = "holm")
#' summary(lm(Y1~time*group, data = dtH0))$coef
#' summary(lm(Y2~time*group, data = dtH0))$coef
#' 
#' res11 <- compMean(Y = cbind(Y1 = dtH0$Y1, Y2 = dtH0$Y2),
#'                   group = dtH0$group, time = dtH0$time, permutation.type = "Y")
#' print(res11)
#' dtWH0 <- melt(dtH0, id.vars = c("group","time"), value.name = "Y")
#' summary(lm(Y~variable*group*time, data = dtWH0))$coef
#'
#' #### 2- H1 under the alternative for the group effect ####
#' dtH1 <- rbind(data.table(Y1 = rnorm(n), Y2 = rnorm(n), group = "C", time = "D1"),
#'               data.table(Y1 = rnorm(n), Y2 = rnorm(n), group = "T", time = "D1"),
#'               data.table(Y1 = rnorm(n), Y2 = rnorm(n), group = "C", time = "W1"),
#'               data.table(Y1 = rnorm(n)+0.5, Y2 = rnorm(n)+0.5, group = "T", time = "W1"))
#'
#' dtH1.mean <- dtH1[,.(meanY1 = mean(Y1)),by=c("group","time")]
#' dtH1.mean[group=="T", diff(meanY1)] - dtH1.mean[group=="C", diff(meanY1)]
#'
#' set.seed(10)
#' res2 <- compMean(Y = cbind(Y1 = dtH1$Y1, Y2 = dtH1$Y2),
#'                  group = dtH1$group, time = dtH1$time, permutation.type = "group")
#' print(res2, method = "max")
#' summary(lm(Y1~time*group, data = dtH1))$coef
#' summary(lm(Y2~time*group, data = dtH1))$coef
#' 
#' res22 <- compMean(Y = cbind(Y1 = dtH1$Y1, Y2 = dtH1$Y2),
#'                   group = dtH1$group, time = dtH1$time, permutation.type = "Y")
#' print(res22, method = "max")
#' dtWH1 <- melt(dtH1, id.vars = c("group","time"), value.name = "Y")
#' summary(lm(Y~variable*group*time, data = dtWH1))$coef
#' 
#' #### 3- H1 under the alternative for the marker effect ####
#' dtH2 <- rbind(data.table(Y1 = rnorm(n), Y2 = rnorm(n), group = "C", time = "D1"),
#'               data.table(Y1 = rnorm(n), Y2 = rnorm(n), group = "T", time = "D1"),
#'               data.table(Y1 = rnorm(n), Y2 = rnorm(n), group = "C", time = "W1"),
#'               data.table(Y1 = rnorm(n) + 0.5, Y2 = rnorm(n) + 1, group = "T", time = "W1"))
#' 
#' res3 <- compMean(Y = cbind(Y1 = dtH2$Y1, Y2 = dtH2$Y2),
#'                  group = dtH2$group, time = dtH2$time, permutation.type = "group")
#' print(res3, method = "max")
#' summary(lm(Y1~time*group, data = dtH2))$coef
#' summary(lm(Y2~time*group, data = dtH2))$coef
#' 
#' res33 <- compMean(Y = cbind(Y1 = dtH2$Y1, Y2 = dtH2$Y2),
#'                   group = dtH2$group, time = dtH2$time, permutation.type = "Y")
#' print(res33)
#' dtWH2 <- melt(dtH2, id.vars = c("group","time"), value.name = "Y")
#' summary(lm(Y~variable*group*time, data = dtWH2))$coef
#'
#' #### 4- Without time ####
#' n <- 500
#' set.seed(10)
#' dtNT <- rbind(data.table(Y1 = rnorm(n), Y2 = rnorm(n), Y3 = rnorm(n), group = "C"),
#'               data.table(Y1 = rnorm(n), Y2 = rnorm(n), Y3 = rnorm(n), group = "T"))
#' res4 <- compMean(Y = cbind(Y1 = dtNT$Y1, Y2 = dtNT$Y2, Y3 = dtNT$Y3),
#'                  group = dtNT$group, permutation.type = "group", n.perm = 1e4)
#' print(res4, method = "max")
#' print(res4, method = "max-step-down")
#' 
#' res44 <- compMean(Y = cbind(Y1 = dtNT$Y1, Y2 = dtNT$Y2),
#'                   group = dtNT$group, permutation.type = "Y")
#' res44
#' 
#' ## comparison to multtest (BiocManager::install("multtest"))
#' library(multtest)
#' mt.maxT(X = t(cbind(Y1 = dtNT$Y1, Y2 = dtNT$Y2, Y3 = dtNT$Y3)),
#'         classlabel = as.numeric(as.factor(dtNT$group))-1)
#' 
#' ## comparison to multcomp
#' library(multcomp)
#' dtNT$group <- as.factor(dtNT$group)
#' lll <- mmm(Y1 = lm(Y1 ~ group, data = dtNT),
#'           Y2 = lm(Y2 ~ group, data = dtNT),
#'           Y3 = lm(Y3 ~ group, data = dtNT))
#' summary(glht(lll, linfct = mlf(mcp(group = "Dunnett"))))
#' 
#' #### 5- With more than two groups ####
#' n <- 500
#' set.seed(10)
#' dtMG <- rbind(data.table(Y1 = rnorm(n), Y2 = rnorm(n), group = "C", time = "D1"),
#'               data.table(Y1 = rnorm(n)+0.1, Y2 = rnorm(n), group = "T1", time = "D1"),
#'               data.table(Y1 = rnorm(n)+0.2, Y2 = rnorm(n), group = "T2", time = "D1"),
#'               data.table(Y1 = rnorm(n), Y2 = rnorm(n), group = "C", time = "W1"),
#'               data.table(Y1 = rnorm(n), Y2 = rnorm(n), group = "T1", time = "W1"),
#'               data.table(Y1 = rnorm(n), Y2 = rnorm(n), group = "T2", time = "W1"))
#' dtMGD1 <- dtMG[time=="D1"]
#' 
#' res5 <- compMean(Y = cbind(Y1 = dtMGD1$Y1, Y2 = dtMGD1$Y2),
#'                  group = dtMGD1$group, permutation.type = "group")
#' print(res5)
#'
#' summary(lm(Y1~group, data = dtMGD1))
#' summary(lm(Y2~group, data = dtMGD1))
#' 
#' res55 <- compMean(Y = cbind(Y1 = dtMG$Y1, Y2 = dtMG$Y2),
#'                   time = dtMG$time, group = dtMG$group, permutation.type = "group")
#' print(res55)
#'
#' summary(lm(Y1~group*time, data = dtMG))
#' summary(lm(Y2~group*time, data = dtMG))
#'
 

## * compMean (code)
##' @export
compMean <- function(Y,
                     group,
                     time = NULL,
                     permutation.type = "group",
                     n.perm = 1e3,
                     na.rm = FALSE,
                     cl = NULL,
                     trace = TRUE){ 
    ## require(gtools)
    requireNamespace("data.table")
    requireNamespace("pbapply")
    n <- NROW(Y)
    p <- NCOL(Y)
    G <- length(unique(group))
    if(!is.matrix(Y)){
        Y <- cbind(Y)
    }
    name <- colnames(Y)
    if(is.null(name)){
        name <- paste0("Y",1:p)
    }
    
    ## ** test args
    permutation.arg <- match.arg(permutation.type, c("group","Y"))
    if(length(group)!=n){
        stop("Argument \'group\' does not has the same length as the number of rows in \'Y\' \n")
    }
    if(!is.null(time)){
        if(length(time)!=n){
            stop("Argument \'time\' does not has the same length as argument \'group\' \n")
        }
        if(length(unique(time))!=2){
            stop("Argument \'time\' must take exactly two different values \n")
        }
    }
    if(is.data.frame(Y)){
        Y <- as.matrix(Y)
    }
    if(!is.numeric(Y)){
        stop("Argument \'Y\' must be numeric \n")
    }
    if(permutation.type == "Y"){
        if(p!=2){
            stop("Argument \'Y\' must should be a matrix with exactly two columns when argument \'permutation.type\' is \"Y\" \n")
        }else{
            Name <- paste0("(",paste0(rev(name), collapse = "-"),")")
        }
        P <- 1
    }else{
        P <- p
        Name <- name
    }
    
    if(!is.null(cl) && trace == FALSE){
        stop("Argument \'cl\' only available when trace is FALSE \n")
    }
    if(any(name == "perm")){
        stop("The \'Y\' argument should not have a column called \'perm\' \n",
             "This name is used internally. \n")
    }

    ## ** contrast function
    Ugroup <- unique(group)
    if(!is.null(time)){
        Utime <- unique(time)
        index.time1 <- which(time==Utime[1])
        index.time2 <- which(time==Utime[2])
        n.time1 <- sum(time==Utime[1])
        n.time2 <- sum(time==Utime[2])
        grid <- cbind(expand.grid(
            strata = Name,
            time = c("T1","T2","T12"),
            group = Ugroup[-1]
        ), estimate = NA, se = NA, statistic = NA)
        grid$collapse <- paste(grid$strata,grid$time,grid$group,sep=".")
        UgroupTime <- expand.grid(Ugroup, Utime)

        calcDiff <- function(Y, group, time, type){

            ls.Y <- lapply(1:NROW(UgroupTime), function(iRow){ ## iRow <- 1
                return(Y[(group==UgroupTime[iRow,1])*(time==UgroupTime[iRow,2])==1,,drop=FALSE])
            })
            names(ls.Y) <- interaction(UgroupTime[,1],UgroupTime[,2])
            
            out <- grid

            for(iG in 2:G){ ## iG <- 2

                    for(iP in 1:P){ ## iP <- 1

                        ## position
                        iPos.T1 <- which(out$collapse == paste0(Name[iP],".T1.",Ugroup[iG]))
                        iPos.T2 <- which(out$collapse == paste0(Name[iP],".T2.",Ugroup[iG]))
                        iPos.T12 <- which(out$collapse == paste0(Name[iP],".T12.",Ugroup[iG]))

                        ## outcome
                        if(permutation.type=="group"){
                            iYref <- list(ls.Y[[paste(Ugroup[1],Utime[1],sep=".")]][,iP],ls.Y[[paste(Ugroup[1],Utime[2],sep=".")]][,iP])
                            iYcomp <- list(ls.Y[[paste(Ugroup[iG],Utime[1],sep=".")]][,iP],ls.Y[[paste(Ugroup[iG],Utime[2],sep=".")]][,iP])
                            if(na.rm){
                                n.iYref <- lapply(iYref,function(x){sum(!is.na(x))})
                                n.iYcomp <- lapply(iYcomp,function(x){sum(!is.na(x))})
                            }else{
                                n.iYref <- lapply(iYref,length)
                                n.iYcomp <- lapply(iYcomp,length)
                            }
                        }else{
                            iYref <- list(ls.Y[[paste(Ugroup[1],Utime[1],sep=".")]],ls.Y[[paste(Ugroup[1],Utime[2],sep=".")]])
                            iYcomp <- list(ls.Y[[paste(Ugroup[iG],Utime[1],sep=".")]],ls.Y[[paste(Ugroup[iG],Utime[2],sep=".")]])
                            if(na.rm){
                                n.iYref <- lapply(iYref,function(x){colSums(!is.na(x))})
                                n.iYcomp <- lapply(iYcomp,function(x){colSums(!is.na(x))})
                            }else{
                                n.iYref <- lapply(iYref,function(x){rep(NROW(x),2)})
                                n.iYcomp <- lapply(iYref,function(x){rep(NROW(x),2)})
                            }
                        }
                        
                        ## Welch t-test
                        if(permutation.type=="group"){
                            out[iPos.T1,"estimate"] <- mean(iYcomp[[1]], na.rm = na.rm) - mean(iYref[[1]], na.rm = na.rm)
                            out[iPos.T1,"se"] <- sqrt(var(iYcomp[[1]], na.rm = na.rm)/n.iYcomp[[1]] + var(iYref[[1]], na.rm = na.rm)/n.iYref[[1]])
                        }else if(permutation.type=="Y"){
                            out[iPos.T1,"estimate"] <- diff(colMeans(iYcomp[[1]]), na.rm = na.rm) - diff(colMeans(iYref[[1]]), na.rm = na.rm)
                            term1.T1 <- var(iYcomp[[1]][,1], na.rm = na.rm)/n.iYcomp[[1]][1] + var(iYcomp[[1]][,2], na.rm = na.rm)/n.iYcomp[[1]][2]
                            term2.T1 <- var(iYref[[1]][,1], na.rm = na.rm)/n.iYref[[1]][1] + var(iYref[[1]][,2], na.rm = na.rm)/n.iYref[[1]][2]
                            out[iPos.T1,"se"] <- sqrt(term1.T1 + term2.T1)
                        }
                        out[iPos.T1,"statistic"] <- out[iPos.T1,"estimate"]/out[iPos.T1,"se"]
                    
                        if(permutation.type=="group"){
                            out[iPos.T2,"estimate"] <- mean(iYcomp[[2]], na.rm = na.rm) - mean(iYref[[2]], na.rm = na.rm)
                            out[iPos.T2,"se"] <- sqrt(var(iYcomp[[2]], na.rm = na.rm)/n.iYcomp[[2]] + var(iYref[[2]], na.rm = na.rm)/n.iYref[[2]])
                        }else if(permutation.type=="Y"){
                            out[iPos.T2,"estimate"] <- diff(colMeans(iYcomp[[2]], na.rm = na.rm)) - diff(colMeans(iYref[[2]], na.rm = na.rm))
                            term1.T2 <- var(iYcomp[[2]][,1], na.rm = na.rm)/n.iYcomp[[2]][1] + var(iYcomp[[2]][,2], na.rm = na.rm)/n.iYcomp[[2]][2]
                            term2.T2 <- var(iYref[[2]][,1], na.rm = na.rm)/n.iYref[[2]][1] + var(iYref[[2]][,2], na.rm = na.rm)/n.iYref[[2]][2]
                            out[iPos.T2,"se"] <- sqrt(term1.T2 + term2.T2)
                        }
                        out[iPos.T2,"statistic"] <- out[iPos.T2,"estimate"]/out[iPos.T2,"se"]

                        out[iPos.T12,"estimate"] <- out[iPos.T2,"estimate"] - out[iPos.T1,"estimate"]
                        out[iPos.T12,"se"] <- sqrt(out[iPos.T1,"se"]^2 + out[iPos.T2,"se"]^2)
                        out[iPos.T12,"statistic"] <- out[iPos.T12,"estimate"]/out[iPos.T12,"se"]
                    }
            }
            return(out)
        }
    }else{
        Utime <- NULL
        grid <- cbind(expand.grid(
            strata = Name,
            group = Ugroup[-1]
        ), estimate = NA, se = NA, statistic = NA)
        grid$collapse <- paste(grid$strata,grid$group,sep=".")
        calcDiff <- function(Y, group, type, ...){
            ls.Y <- lapply(Ugroup, function(iG){
                if(permutation.type == "group"){
                    return(Y[(group==iG),,drop=FALSE])
                }else if(permutation.type == "Y"){
                    return(Y[(group==iG),2,drop=FALSE]-Y[(group==iG),1,drop=FALSE])
                }
            })
            names(ls.Y) <- Ugroup
            
            out <- grid
            for(iG in 2:G){

                for(iP in 1:P){

                    ## position
                    iPos <- which(out$collapse == paste0(Name[iP],".",Ugroup[iG]))

                    ## outcome
                    iYref <- ls.Y[[Ugroup[1]]][,iP]
                    iYcomp <- ls.Y[[Ugroup[iG]]][,iP]

                    ## Welch t-test
                    out[iPos,"estimate"] <- mean(iYcomp) - mean(iYref)
                    out[iPos,"se"] <- sqrt(var(iYcomp)/length(iYcomp) + var(iYref)/length(iYref))
                    out[iPos,"statistic"] <- out[iPos,"estimate"]/out[iPos,"se"]
                    
                }
            }
            return(out)
        }
    }
    
    ## ** compute means
    out <- list(estimate = calcDiff(Y = Y, group = group, time = time, type = permutation.type),
                time = Utime,
                Yname = name,
                group = Ugroup,
                permutation.type = permutation.type,
                n.perm = n.perm)
    ## class(out$estimate)

    ## ** run permutation

    fct.apply <- switch(as.character(trace),
                        "TRUE" = pbapply::pblapply,
                        "FALSE" = lapply)
    ls.out <- do.call(fct.apply, args = list(1:n.perm, function(iPerm){## iPerm <- 1
        if(permutation.type == "group"){
            if(!is.null(time)){
                newgroup <- rep(NA, length = n)
                newgroup[index.time1] <- sample(group[time==Utime[1]], size = n.time1, replace = FALSE)
                newgroup[index.time2] <- sample(group[time==Utime[2]], size = n.time2, replace = FALSE)
            }else{
                newgroup <- sample(group, size = n, replace = FALSE)
            }
            iOut <- calcDiff(Y = Y, group = newgroup, time = time)
            return(cbind(out$estimate, perm = iPerm, perm.estimate = iOut$estimate, perm.se = iOut$se, perm.statistic = iOut$statistic))
        }else if(permutation.type == "Y"){
            newY <- Y
            indexExchange <- rbinom(n, size = 1, prob = 0.5)
            newY[indexExchange==1,1] <- Y[indexExchange==1,2]
            newY[indexExchange==1,2] <- Y[indexExchange==1,1]
            iOut <- calcDiff(Y = newY, group = group, time = time, type = permutation.type)
            
            return(cbind(out$estimate, perm = iPerm, perm.estimate = iOut$estimate, perm.se = iOut$se, perm.statistic = iOut$statistic))
        }
    }, cl = cl))
    
    out$permutation <- as.data.table(do.call(rbind, ls.out))
    ## p-value
    if(is.null(time)){out$permutation[,c("time") := "T1"]}
    out$p.value <- out$permutation[,.(estimate = mean(abs(perm.estimate)>abs(estimate)),
                                      statistic = mean(abs(perm.statistic)>abs(statistic))),
                                   by = c("strata","group","time")]
    ## adjusted p-value (holm)
    out$p.value$adjHolm.statistic <- p.adjust(out$p.value$statistic, method = "holm")

    ## adjusted p-value (max)
    iStatistic.max <- abs(out$estimate$statistic)
    iRank <- rank(iStatistic.max)  ## find the rank of the observed statistic 

    permutationW <- out$permutation[,setNames(as.list(abs(perm.statistic)),collapse), by = "perm"] ## reshape the permutation statistics
    name.test <- setdiff(names(permutationW),"perm")
    permutation.max <- apply(permutationW[,.SD,.SDcols = name.test],1,max)
    
    out$p.value$adjMax.statistic <- NA
    out$p.value$adjMaxDown.statistic <- NA
    p.test <- switch(permutation.type,
                     "group" = length(unique(out$p.value$time))*p*(G-1),
                     "Y" = length(unique(out$p.value$time)))
    for(iTest in 1:p.test){ ## iTest <- 1
        iIndex <- which(iRank == (p.test - iTest + 1) ) ## will not work with ex-aequo

        ## always take the largest
        out$p.value$adjMax.statistic[iIndex] <- max(c(mean(permutation.max > abs(iStatistic.max[iIndex])),
                                                      out$p.value$adjMax.statistic),
                                                    na.rm = TRUE)
        ## take the iTest-th largest
        iRemaining <- name.test[which(is.na(out$p.value$adjMaxDown.statistic))]
        iPermutation.max <- apply(permutationW[,.SD,.SDcols = iRemaining],1,max)
        out$p.value$adjMaxDown.statistic[iIndex] <- max(c(mean(iPermutation.max > abs(iStatistic.max[iIndex])),
                                                          out$p.value$adjMaxDown.statistic),
                                                        na.rm = TRUE)

    }
    ## output
    class(out) <- "permTest"
    return(out)
}

## * print.permTest
##' @export
print.permTest <- function(x, method = "max-step-down", setkeyv = NULL, ...){
    time <- x$time
    ref <- x$group[1]
    other <- setdiff(x$group,ref)
    if(length(other)>1){
        label.other <- "."
    }else{
        label.other <- other[1]
    }
    
    if(x$permutation.type=="group"){
        
        cat("\n     Comparison between ",length(x$group)," groups with respect to each marker Y\n\n",sep="")
        cat(" - groups:  ",paste(other,collapse="/")," (alternative) vs. ",ref," (reference) \n",sep="")

        if(!is.null(time)){
            Y.T1ref <- paste0("Y(",x$time[1],",",ref,")")
            Y.T2ref <- paste0("Y(",x$time[2],",",ref,")")
            Y.T1comp <- paste0("Y(",x$time[1],",",label.other,")")
            Y.T2comp <- paste0("Y(",x$time[2],",",label.other,")")

            cat(" - times :\n",
                "     at the first time  (",x$time[1],")   :  ",Y.T1comp," - ",Y.T1ref,"\n",
                "     at the second time (",x$time[2],")   :  ",Y.T2comp," - ",Y.T2ref,"\n",
                "     temporal evolution (",x$time[2],"-",x$time[1],"): (",Y.T2comp," - ",Y.T2ref,") - (",Y.T1comp," - ",Y.T1ref,") \n",
                sep="")
        }else{
            Y.G1 <- paste0("Y(",ref,")")
            Y.G2 <- paste0("Y(",label.other,")")
            cat(" - test  : ",Y.G2," - ",Y.G1,"\n",sep="")
        }
        
    }else if(x$permutation.type=="Y"){
    
        cat("\n     Comparison of the group effect between two markers Y1 and Y2\n\n")
        cat(" - groups:  ",paste(other,collapse="/")," (alternative) vs. ",ref," (reference) \n",sep="")

        if(!is.null(time)){
            Y1.T1ref <- paste0("Y1(",x$time[1],",",ref,")")
            Y1.T2ref <- paste0("Y1(",x$time[2],",",ref,")")
            Y1.T1comp <- paste0("Y1(",x$time[1],",",label.other,")")
            Y1.T2comp <- paste0("Y1(",x$time[2],",",label.other,")")
            Y2.T1ref <- paste0("Y2(",x$time[1],",",ref,")")
            Y2.T2ref <- paste0("Y2(",x$time[2],",",ref,")")
            Y2.T1comp <- paste0("Y2(",x$time[1],",",label.other,")")
            Y2.T2comp <- paste0("Y2(",x$time[2],",",label.other,")")

            cat(" - times :\n",
                "     at the first time  (",x$time[1],")   :  (",Y2.T1comp," - ",Y2.T1ref,") - (",Y1.T1comp," - ",Y1.T1ref,")\n",
                "     at the second time (",x$time[2],")   :  (",Y2.T2comp," - ",Y2.T2ref,") - (",Y1.T2comp," - ",Y1.T2ref,")\n",
                "     temporal evolution (",x$time[2],"-",x$time[1],"): [(",Y2.T2comp," - ",Y2.T2ref,") - (",Y1.T2comp," - ",Y1.T2ref,")] \n",
                "                                  - [(",Y2.T1comp," - ",Y2.T1ref,") - (",Y1.T1comp," - ",Y1.T1ref,")] \n",
                sep="")

        }else{
            Y1.ref <- paste0("Y1(",ref,")")
            Y2.ref <- paste0("Y2(",ref,")")
            Y1.comp <- paste0("Y1(",label.other,")")
            Y2.comp <- paste0("Y2(",label.other,")")
            cat(" - test  : (",Y2.comp," - ",Y2.ref,") - (",Y1.comp," - ",Y1.ref,")\n",sep ="")

        }
    }
    
    cat("\n")

    if(!is.null(time)){
        table <- cbind(x$estimate[,c("strata","group","time","estimate","se")],
                       p.value = x$p.value$statistic)
        names(table)[names(table)=="group"] <- "alternative"
    }else{
        table <- cbind(x$estimate[,c("strata","group","estimate","se")],
                       p.value = x$p.value$statistic)
        names(table)[names(table)=="group"] <- "alternative"
    }
    if(method == "holm"){
        table$adj.p.value <- x$p.value$adjHolm.statistic
    }else if(method == "max"){
        table$adj.p.value <- x$p.value$adjMax.statistic
    }else if(method == "max-step-down"){
        table$adj.p.value <- x$p.value$adjMaxDown.statistic
    }else stop("Incorrect \'method\' argument, must be \"holm\", \"max\", or \"max-step-down\". \n")
    
    table <- as.data.table(table)
    if(!is.null(time)){
        setkeyv(table, c("strata","time"))
        table$time <- factor(table$time, levels = c("T1","T2","T12"), labels = c(x$time[1],x$time[2],paste0(x$time[2],"-",x$time[1])))
    }else{
        setkeyv(table, c("strata"))
    }
    if(x$permutation.type=="group"){
        table$strata <- as.character(table$strata)
        table$strata[duplicated(table$strata)] <- ""
        names(table)[1] <- "Y"
    }else if(x$permutation.type == "Y"){
        table$strata <- NULL
        table$adj.p.value <- NULL
    }

    ## display
    if(!is.null(setkeyv)){
        data.table::setkeyv(table, setkeyv)
    }
    print(table, row.names = FALSE)
    cat("(based on ",x$n.perm," repetitions)\n",sep="")
    if(x$permutation.type=="group"){
        cat("(adjustment for multiple comparisons: ",method,")\n",sep="")
    }
    
    ## export
    return(invisible(table))
}


######################################################################
### compMean.R ends here
