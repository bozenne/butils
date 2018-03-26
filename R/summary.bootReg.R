### summary.bootReg.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 21 2017 (17:49) 
## Version: 
## Last-Updated: mar 26 2018 (14:21) 
##           By: Brice Ozenne
##     Update #: 296
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * summary.bootReg
#' @title Display the result fo the bootstrap computation
#' @description Display the result fo the bootstrap computation
#' 
#' @param object object obtained with the function \code{bootReg}
#' @param p.value should the p.value be computed? Can be time consuming.
#' @param conf the confidence level of the confidence intervals.
#' @param type the method used to compute the confidence intervals.
#' See the documentation of \code{boot::boot.ci}.
#' @param index for which parameters confidence intervals/p.values should be computed? Default is for all parameters.
#' @param print should the summary be displayed in the terminal?
#' @param n.subset if not \code{NULL} display in addition intervals/p.values for subset of the bootstrap samples.
#' @param ... not used
#'
#' @details This function calls the \code{boot.ci} from the boot package
#' to compute the confidence intervals.
#'
#' P.value are computed by finding the confidence level
#' for which one bound of the confidence interval touches 0.
#'
#' When the number of bootstrap sample is too low, the function (in fact \code{boot.ci}) will return warnings like
#' "extreme order statistics used as endpoints".
#' 
#' @method summary bootReg
#' @export
summary.bootReg <- function(object, p.value = TRUE,
                            conf = 0.95, type = NULL, index = NULL,
                            print = TRUE, n.subset = NULL, ...){

    if(is.null(index)){
        index <- 1:length(object$estimate)
    }else if(is.character(index)){
        index <- match(index,object$estimate)
    }

    if(is.null(type)){
        if(!is.null(object$FCT.stdError)){
            type <- "stud"
        }else if(!is.null(object$FCT.iid)){
            type <- "bca"
        }else{
            type <- "perc"
        }
    }
    
    ### ** Convert to boot
    boot.object <- as.boot(object, index = index)

    ### ** Compute confidence intervals and p.values
    resStat.full <- .calcCIP(boot.object, conf = conf, type = type, index = 1:length(index),
                             p.value = p.value,
                             stdError = object$stdError[index],
                             iid = object$iid[,index,drop=FALSE],
                             boot.stdError = object$boot.stdError[,index,drop=FALSE])
    
    ### ** Display results
    rowM <- cbind(estimate = object$estimate[index],
                  estimate.boot = apply(object$boot.estimate[,index,drop=FALSE],2,median,na.rm = TRUE),
                  resStat.full,
                  nBoot.effective = colSums(!is.na(object$boot.estimate[,index,drop=FALSE])))

    if(print==TRUE){
        print(object$call)
        cat("\n")
        print(rowM)
        cat("\n")
        if(is.null(object$FUN.resample)){
            cat("Non-parametric bootstrap with ",as.double(object$n.boot)," replications \n", sep = "")
        }else if(object$FUN.resample %in% c("simulate","sim")){
            cat("Parametric bootstrap with ",as.double(object$n.boot)," replications \n", sep = "")            
        } else{
            cat("bootstrap with ",as.double(object$n.boot)," replications \n", sep = "")
        }

        label.type <- switch(type,
                             "norm" = "normal",
                             "basic" = "basic",
                             "stud" = "studentized",
                             "perc" = "percentile",
                             "bca" = "BCa")
        cat("Confidence intervals estimated using ",label.type," boostrap - see help(boot.ci) \n", sep = "")        
    }

### ** By subset
    if(!is.null(n.subset)){ # seq.length.out <- 5
        cat("\n\n") 
        
        seqIndex <- round(seq(1,boot.object$R+1, length.out = n.subset+1))
        M.index <- cbind(seqIndex[-length(seqIndex)],(seqIndex-1)[-1])
        n.Index <- NROW(M.index)
  
        resSubset <- lapply(1:n.Index, function(iSub){ # iSub <- 1
            iIndex <- M.index[iSub,1]:M.index[iSub,2]
                
            iBoot.object <- boot.object
            iBoot.object$t <- iBoot.object$t[iIndex,,drop=FALSE]
            iBoot.object$strata <- iBoot.object$strata[iIndex]
            iBoot.object$weights <- iBoot.object$weight[iIndex]
            iBoot.object$R <- NROW(iBoot.object$t)
            resStat.full <- .calcCIP(iBoot.object, conf = conf, type = type, index = index,
                                     stdError = object$stdError,
                                     iid = object$iid[iIndex,,drop=FALSE],
                                     boot.stdError = object$boot.stdError[iIndex,,drop=FALSE])
            return(resStat.full)
        })
            
        print(resSubset)
        cat("for subsets \n")
    }

### ** export
    return(invisible(rowM))
}

## * .calcCIP
.calcCIP <- function(boot.object, conf, type, index,
                     p.value,
                     stdError, boot.stdError, iid,
                     ...){

    type <- match.arg(type, c("norm","basic","stud","perc","bca"))
    slot.boot.ci <- switch(type,
                           "norm" = "normal",
                           "basic" = "basic",
                           "stud" = "student",
                           "perc" = "percent",
                           "bca" = "bca")
    
    n.boot <- boot.object$R
    grid <- seq(0,by=1/n.boot,length.out=n.boot)
    name.estimate <- names(boot.object$t0)
    n.estimate <- length(name.estimate)

    index.lowerCI <- switch(type,
                            "norm" = 2,
                            "basic" = 4,
                            "stud" = 4,
                            "perc" = 4,
                            "bca" = 4)
    index.upperCI <- switch(type,
                            "norm" = 3,
                            "basic" = 5,
                            "stud" = 5,
                            "perc" = 5,
                            "bca" = 5)

    ### ** functions
    calcCI <- function(boot.object, conf, type, index,
                       stdError, boot.stdError, iid,
                       slot.boot.ci, name.estimate){
        ls.out <- lapply(index, function(iP){ # iP <- 1
          
            if(type == "norm"){
                resBoot.ci <- boot::boot.ci(boot.object, conf = conf, type = type, index = iP)
            }else{
                iBoot.object <- boot.object
                iBoot.object$t0 <- c(boot.object$t0[iP],stdError[iP])
                iBoot.object$t <- cbind(boot.object$t[,iP,drop=FALSE],boot.stdError[,iP,drop=FALSE])

                resBoot.ci <- boot::boot.ci(iBoot.object, conf = conf, type = type,
                                            var.t0 = stdError, var.t = boot.stdError[,iP],
                                            t0 = boot.object$t0[iP], t = boot.object$t[,iP],
                                            L = iid[,iP])
            }
            out <- resBoot.ci[[slot.boot.ci]][index.lowerCI:index.upperCI]
            return(setNames(out,c("lower","upper")))
        })
        out <- do.call(rbind,ls.out)
        rownames(out) <- name.estimate[index]
        return(out)
    }

    
        
    calcPvalue <- function(boot.object, conf, type, index,
                           stdError, boot.stdError, iid,
                           slot.boot.ci, name.estimate){
       
        out <- sapply(index, function(iP){ # iP <- 1
            if(abs(boot.object$t0[iP])<.Machine$double.eps^0.5){
                return(0)                
            }else{               
                iBoot.object <- boot.object
                iBoot.object$t0 <- c(boot.object$t0[iP],stdError[iP])
                iBoot.object$t <- cbind(boot.object$t[,iP,drop=FALSE],boot.stdError[,iP,drop=FALSE])

                ## search confidence level such that quantile of CI which is close to 0
                p.value <- boot2pvalue(x = boot.object$t[,iP],
                                       estimate = boot.object$t0[iP],
                                       alternative = "two.sided",
                                       FUN.ci = function(x, p.value, sign.estimate, ...){                                           
                                           side.CI <- c(index.lowerCI,index.upperCI)[2-sign.estimate]
                                           if(type == "norm"){
                                               iRes <- boot::boot.ci(iBoot.object, conf = 1-p.value, type = type,
                                                                     index = 1:2)[[slot.boot.ci]][side.CI]                            
                                           }else{
                                               iRes <- boot::boot.ci(iBoot.object, conf = 1-p.value, type = type, index = 1:2,
                                                                     var.t0 = stdError[iP], var.t = boot.stdError[,iP],
                                                                     t0 = boot.object$t0[iP], t = boot.object$t[,iP],
                                                                     L = iid[,iP])[[slot.boot.ci]][side.CI]                           
                                           }
                                       })
                ##boot::boot.ci(iBoot.object, conf = 0.89, type = type, index = 1:2)
                return(p.value)
            }
            out <- setNames(out, name.estimate[index])
            return(out)
        })
    }
    
    ### ** compute p values and CI for each subset
    resCI <- calcCI(boot.object = boot.object, conf = conf, type = type, index = index,
                    stdError = stdError, boot.stdError = boot.stdError, iid = iid,
                    slot.boot.ci = slot.boot.ci, name.estimate = name.estimate)

    if(p.value){
        resP <- calcPvalue(boot.object = boot.object,
                           conf = conf,
                           type = type,
                           index = index,
                           stdError = stdError,
                           boot.stdError = boot.stdError,
                           iid = iid,
                           slot.boot.ci = slot.boot.ci,
                           name.estimate = name.estimate)
        resCI <- cbind(resCI,p.value = resP)
    }

### ** export
    return(resCI)
}



##----------------------------------------------------------------------
### summary.bootReg.R ends here
