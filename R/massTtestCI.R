### massTtestCI.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 29 2018 (09:33) 
## Version: 
## Last-Updated: dec  6 2018 (17:29) 
##           By: Brice Ozenne
##     Update #: 57
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


## * massTtestCI (documentation)
##' @title Mass t-tests with conditional inference
##' @description Perform t-tests for all variables or only for variables showing a sufficiently large effect.
##' In the latter case, conditional confidence intervals are used instead of the "normal" confidence intervals
##' to account for the selection step. 
##' @name massTtestCI
##'
##' @param X [numeric matrix] the design matrix for the individuals of the control group. 
##' @param Y [numeric matrix] the design matrix for the individuals of the experimental group.
##' @param threshold [numeric] the selection threshold for the difference in mean between the groups.
##' @param conf.level [numeric, 0-1] Confidence level.
##' @param method.p.adjust [character] Method used to adjust for multiple comparisons.
##' See argument method in \code{p.adjust} for more details.
##' @param ... additional arguments passed to the function \code{conditionalCI}.
##' 
##' @return A data.frame where each site corresponds to a row and with the following columns: \itemize{
##' \item site: the index of the site
##' \item estimate: the estimated mean difference
##' \item sigma: the estimated standard error of the estimated mean difference
##' \item statistic: the value of the test statistic
##' \item p.value: the (unadjusted) p-value 
##' \item adjusted.p.value: the p-value adjusted for multiple comparisons over all sites
##' \item select: whether the site has been selected (i.e. signal above the threshold)
##' \item lower.CI: lower bound of the post-selection confidence interval (unadjusted)
##' \item upper.CI: upper bound of the post-selection confidence interval (unadjusted)
##' \item adjusted.lower.CI: lower bound of the post-selection confidence interval (adjusted for multiple comparisons over all selected sites)
##' \item adjusted.upper.CI: upper bound of the post-selection confidence interval (adjusted for multiple comparisons over all selected sites)
##' }

## * massTtestCI (examples)
##' @rdname massTtestCI
##' @examples
##'
##' if(require(sp) & require(gstat)){
##'
##' ## settings
##' n.samples <- 20
##' n.sites <- 30*30
##'
##' coords <- expand.grid(x = 1:sqrt(n.sites),
##'                      y = 1:sqrt(n.sites))
##' coords$mu <- 0
##' coords$mu[coords$x<=5 & coords$y<=5] <- 1
##'
##' vario <- vgm(psill = 1, range = 2, model='Exp')
##' 
##' ## simulate data
##' set.seed(10)
##' 
##' X <- sim2Dimage(n = n.samples,
##'                 coords = coords[,c("x","y")],
##'                 mu = 0,
##'                 vgm = vario)$X
##' 
##' Y <- sim2Dimage(n = n.samples,
##'                 coords = coords[,c("x","y")],
##'                 mu = coords$mu,
##'                 vgm = vario)$X
##'
##' coords$meanX <- colMeans(X)
##' coords$meanY <- colMeans(Y)
##'
##' spdf <- coords
##' spdf$effect <- coords$meanY - coords$meanX
##' gridded(spdf) <- ~x+y
##' effect.range <- range(spdf[["effect"]])
##' 
##' plot(spdf[,"effect"], zlim = effect.range)
##'
##' res <- massTtestCI(X, Y, threshold = 0.5, method.CI = "shortest2")
##'
##' ## traditional analysis
##' spdf$p.value <- res$p.value
##' spdf$adjusted.p.value <- res$adjusted.p.value
##' 
##' plot(spdf[,"effect"], zlim = effect.range,
##' main = "traditional t-test")
##' points(spdf[spdf$adjusted.p.value<0.05,"p.value"], col = "green", lwd = 3)
##'
##' ## post selection CI
##' spdf$selected <- res$selected
##' spdf$lower.CI <- res$lower.CI
##' spdf$adjusted.lower.CI <- res$adjusted.lower.CI
##' 
##' plot(spdf[spdf$selected>0,"lower.CI"], zlim = effect.range,
##'      main = "post-selection unadjusted")
##' points(spdf[which(spdf$lower.CI>0),"lower.CI"], col = "green", lwd = 3)
##'
##' plot(spdf[spdf$selected>0,"adjusted.lower.CI"], zlim = effect.range,
##'      main = "post-selection adjusted")
##' points(spdf[which(spdf$adjusted.lower.CI>0),"lower.CI"], col = "green", lwd = 3)
##' }
##' 

## * massTtestCI (code)
##' @rdname massTtestCI
##' @export
massTtestCI <- function(X, Y, threshold, conf.level = 0.95, method.p.adjust = "bonferroni",
                        ...){

    alpha <- 1 - conf.level
    n.sites <- NCOL(X)
    if(NCOL(Y)!=n.sites){
        stop("Arguments \'X\' and \'Y\' must have the same number of columns \n")
    }
    out <- data.frame(site = 1:n.sites)
    
    ## ** mass t-tests
    ls.ttest <- lapply(1:n.sites, function(i){stats::t.test(X[,i],Y[,i])})

    
    out$estimate <- sapply(ls.ttest,function(x){as.double(diff(x[["estimate"]]))})
    out$statistic <- sapply(ls.ttest,function(x){as.double(x[["statistic"]])})
    out$sigma <- abs(out$estimate)/abs(out$statistic)
    out$p.value <- sapply(ls.ttest,"[[","p.value")
    
    ## ** adjusted p-values
    out$adjusted.p.value <- stats::p.adjust(out$p.value, method = method.p.adjust)

    ## ** conditional inference
    ## not adjusted for multiple comparisons
    iCCI <- try(conditionalCI(theta = out$estimate,
                              sigma = out$sigma,
                              threshold = threshold,
                              conf.level = conf.level,
                              trace = FALSE,
                              ...))
        
    if(!inherits(iCCI,"try-error")){
        out$selected <- !is.na(stats::confint(iCCI)$value)
        out$lower.CI <- iCCI$CI[,"lower"]
        out$upper.CI <- iCCI$CI[,"upper"]
        out$significant.CI <- FALSE
        out$significant.CI[which(iCCI$CI[,"lower"]>0)] <- TRUE
        out$significant.CI[which(iCCI$CI[,"upper"]<0)] <- TRUE

        n.selected <- sum(out$selected)
        ## not adjusted for multiple comparisons using bonferroni
        adjusted.iCCI <- try(conditionalCI(theta = out$estimate,
                                           sigma = out$sigma,
                                           threshold = threshold,
                                           conf.level = 1 - alpha/n.selected,
                                           trace = FALSE,
                                           ...))
        if(!inherits(adjusted.iCCI,"try-error")){
            out$adjusted.lower.CI <- adjusted.iCCI$CI[,"lower"]
            out$adjusted.upper.CI <- adjusted.iCCI$CI[,"upper"]
            out$adjusted.significant.CI <- FALSE
            out$adjusted.significant.CI[which(adjusted.iCCI$CI[,"lower"]>0)] <- TRUE
            out$adjusted.significant.CI[which(adjusted.iCCI$CI[,"upper"]<0)] <- TRUE
        }else{
            out$adjusted.lower.CI <- NA
            out$adjusted.upper.CI <- NA
            out$adjusted.significant.CI <- NA
        }

    }else{
        out$lower.CI <- NA
        out$upper.CI <- NA
        out$significant.CI <- NA
        out$adjusted.lower.CI <- NA
        out$adjusted.upper.CI <- NA
        out$adjusted.significant.CI <- NA
    }
    

    ## ** export
    return(out)    
}

######################################################################
### massTtestCI.R ends here
