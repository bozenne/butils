### lvm2cor.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 17 2021 (20:28) 
## Version: 
## Last-Updated: sep  7 2021 (11:47) 
##           By: Brice Ozenne
##     Update #: 124
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * lvmCov2Cor (documentation)
##' @name lvm2cor
##' @title Covariance To Correlation Coefficient in LVM
##' @description Convert a covariance coefficient into a correlation coefficient in a lvm.
##' Uses a delta method to compute the variance, confidence interval, and p-value associated with the test of no correlation.
##' @param object a latent variable model (\code{lvmfit} object).
##' @param var1 the name of a latent or endogenous variable.
##' @param var2 the name of another latent or endogenous variable.
##' @param level the confidence level used for the confidence interval.
##' @param cluster the grouping variable relative to which the observations are iid.
##' @param ssc should the standard errors of the coefficients be corrected for small sample bias?
##'
##' @return A data.frame containing in the fourth row the estimated correlation coefficient, standard error, confidence intervals, and p.value.
##' 
##' @examples
##' 
##' #### 0 - Simulate some data ####
##' library(lava)
##' mSim <-lvm(c(PEQ_poslife,PEQ_posself,PEQ_posmood,PEQ_possoc,PEQ_posbehav)~lv.peq,
##'            c(MEQ_mystical,MEQ_mood) ~ 1*lv.meq,
##'            c(MEQ_timespace,MEQ_ineffability) ~ lv.meq,
##' lv.meq[0:2]~1,
##' lv.peq[0:0.25]~1)
##' covariance(mSim) <- lv.peq ~ lv.meq
##' covariance(mSim) <- MEQ_timespace~MEQ_ineffability
##' latent(mSim) <- ~lv.peq+lv.meq
##'
##' set.seed(10)
##' GS <- sim(mSim, 1e5, latent = TRUE) ## not observed
##' set.seed(10)
##' d <- sim(mSim, 1e3, latent = FALSE) ##  observed
##'
##' #### Example 1: Fit a LVM with a single latent variable ####
##' m1 <-lvm(c(PEQ_poslife,PEQ_posself,PEQ_posmood,PEQ_possoc,PEQ_posbehav)~lv.peq)
##' latent(m1) <- ~lv.peq
##' e1 <- estimate(m1, d)
##' 
##' ### 1.1 Correlation between endogenous
##' ## lava output
##' keep.col <- c("PEQ_poslife","PEQ_posself")
##' attr(predict(e1),"cond.var")[keep.col,keep.col]
##' cov2cor(attr(predict(e1),"cond.var")[keep.col,keep.col])
##' 
##' ## proposed method
##' lvmCov2Cor(e1, var1 = "PEQ_poslife", var2 = "PEQ_posself")
##'
##' #### Example 2: Fit a LVM with two latent variables ####
##' m2 <- lvm(c(PEQ_poslife,PEQ_posself,PEQ_posmood,PEQ_possoc,PEQ_posbehav)~lv.peq,
##'          c(MEQ_mystical,MEQ_mood) ~ 1*lv.meq,
##'          c(MEQ_timespace,MEQ_ineffability) ~ lv.meq)
##' covariance(m2) <- lv.peq ~ lv.meq
##' covariance(m2) <- MEQ_timespace~MEQ_ineffability
##' latent(m2) <- ~lv.peq + lv.meq
##' e2 <- estimate(m2, d)
##'
##' ### 2.1: Correlation between endogenous
##' ## approximated true value
##' c("var.meq" = var(GS$MEQ_timespace), "var.peq" = var(GS$MEQ_ineffability),
##'   "cov" = cov(GS$MEQ_timespace,GS$MEQ_ineffability),
##'   "cor" = cor(GS$MEQ_timespace,GS$MEQ_ineffability))
##'
##' ## estimate value
##' coef(e2, type = 9)["MEQ_timespace~~MEQ_ineffability",]
##' 
##' ## proposed method
##' lvmCov2Cor(e2, var1 = "MEQ_timespace", var2 = "MEQ_ineffability")
##' lvmCov2Cor(e2, var1 = "MEQ_timespace", var2 = "MEQ_ineffability",
##'            robust = TRUE)
##' 
##' ## using lava 
##' estimate(e2, f = function(x){
##' 
##' sigma1 <- x["MEQ_timespace~~MEQ_timespace"]
##' sigma2 <- x["MEQ_ineffability~~MEQ_ineffability"]
##' sigma12 <- x["MEQ_timespace~~MEQ_ineffability"]
##' lambda1 <- x["MEQ_timespace~lv.meq"]
##' lambda2 <- x["MEQ_ineffability~lv.meq"]
##' tau <- x["lv.meq~~lv.meq"]
##' 
##' Sigma1 <- sigma1 + lambda1^2*tau
##' Sigma2  <- sigma2 + lambda2^2*tau
##' 
##' return(c(sigma12/sqrt(Sigma1*Sigma2),
##'          (sigma12+lambda1*lambda2*tau)/sqrt(Sigma1*Sigma2)))
##' })
##' 
##' ### 2.1: Correlation between latent variables
##' ## approximated true value
##' c("var.meq" = var(GS$lv.meq), "var.peq" = var(GS$lv.peq),
##'   "cov" = cov(GS$lv.meq,GS$lv.peq), "cor" = cor(GS$lv.meq,GS$lv.peq))
##' 
##' ## estimate value
##' coef(e2, type = 9)["lv.peq~~lv.meq",]
##' 
##' ## Delta method via lvmCov2Cor
##' lvmCov2Cor(e2, var1 = "lv.meq", var2 = "lv.peq")
##' lvmCov2Cor(e2, var1 = "lv.meq", var2 = "lv.peq", robust = TRUE)
##' 
##' ## using lava
##' estimate(e2, function(x){
##' a <- x["lv.meq~~lv.meq"]
##' b <- x["lv.peq~~lv.peq"]
##' c <- x["lv.peq~~lv.meq"]
##'   c(var.meq = a,
##'    var.peq = b,
##'    cov = c,
##'    cor = c/sqrt(a*b))
##' })
##'
##' ## not sure what is
##' keep.col <- c("lv.meq","lv.meq")
##' attr(predict(e2, x = manifest(e2), y = latent(e2)),"cond.var")
##' cov2cor(attr(predict(e2, x = manifest(e2), y = latent(e2)),"cond.var"))
##' 
##' 


## * lvmCov2Cor (code)
##' @rdname lvm2cor
##' @export
lvmCov2Cor <- function(object, var1, var2, null = 0, level = 0.95, robust = FALSE, ssc = FALSE, cluster = NULL){
    requireNamespace(numDeriv)

    ## ** extract parameters and their variance-covariance structure from objects
    object.param <- stats::coef(object)
    name.param <- names(object.param)
    
    if(robust){
        if(ssc){
            requireNamespace("lavaSearch2")
            if(!is.null(cluster)){
                object.iid <- lavaSearch2::iid2(object, cluster = cluster)
            }else{
                object.iid <- lavaSearch2::iid2(object)
            }
        }else{
            requireNamespace("lava")
            if(!is.null(cluster)){
                object.iid <- lava::iid(object, id = cluster)
            }else{
                object.iid <- lava::iid(object)
            }
        }
        object.vcov <- crossprod(object.iid)
    }else{
        if(!is.null(cluster)){
            stop("Arugment \'robust\' must be true when argument \'cluster\' is not NULL. \n")
        }
        if(ssc){
            object.vcov <- lavaSearch2::vcov2(object)
        }else{
            object.vcov <- stats::vcov(object)
        }
    }

    ## ** check arguments
    alpha <- 1-level
    if(any("cor" %in% name.param)){
        stop("No coefficient should be called \"cor\". This name is used internally. \n")
    }
    if(any("cov" %in% name.param)){
        stop("No coefficient should be called \"cov\". This name is used internally. \n")
    }
    if(any("var1" %in% name.param)){
        stop("No coefficient should be called \"var1\". This name is used internally. \n")
    }
    if(any("var2" %in% name.param)){
        stop("No coefficient should be called \"var2\". This name is used internally. \n")
    }
    if(var1 %in% lava::vars(object) == FALSE){
        stop("Argument \'var1\' does not correspond to a variable of the LVM. \n")
    }
    if(var2 %in% lava::vars(object) == FALSE){
        stop("Argument \'var2\' does not correspond to a variable of the LVM. \n")
    }
    if(var1 == var2){
        stop("Argument \'var2\' incorrect: should not be the same name as argument \'var1\'. \n")
    }

    ## ** prepare function that outputs correlation
    ## measurement model: Y = \alpha + \Gamma X + \Lambda \eta + \varepsilon
    ## structural model: \eta = \nu + K X + B \eta + \xi
    ## overall model: Z = a + b X + c Z + \nu
    ## so Z = (1-c)^{-1}(a+b) + (1-c)^{-1} \nu
    
    FUN.cor <- function(p){ ## p <- coef(object)
        iMoment <- lava::moments(object, p = p)
        iIndex <- lava::index(object)

        ## residual variance for latent and endogenous variables
        Sigma.nu <- iMoment$P
        Sigma.nu[iIndex$exo.idx, iIndex$exo.idx] <- 0  ## condition on exogenous variables
        ## link function with parents
        M.c <- iMoment$A 
        M.1mcM1 <- solve(diag(1, nrow = NROW(M.c), ncol = NCOL(M.c))-t(M.c)) ## inverse link function

        if(any(abs(M.1mcM1 - iMoment$IAi)>1e-10)){
            warning("Something went wrong when inverting the link matrix. \n",
                     "Do not match lava results")
        }

        ## variance conditional on X = (1-c)^{-1} \nu t((1-c)^{-1} )
        Sigma.conditional <- M.1mcM1 %*% Sigma.nu %*% t(M.1mcM1)
        sigma.var1 <- Sigma.conditional[var1,var1]
        sigma.var2 <- Sigma.conditional[var2,var2]

        ## direct covariance
        sigma.var12_direct <- Sigma.nu[var1,var2]

        ## total covariance
        sigma.var12_total <- Sigma.conditional[var1,var2]
        
        ## result
        out <- c(var1 = sigma.var1,
                 var2 = sigma.var2,
                 cov12_direct = sigma.var12_direct,
                 cov12_total = sigma.var12_total,
                 cor12_direct = sigma.var12_direct/sqrt(sigma.var1*sigma.var2),
                 cor12_total = sigma.var12_total/sqrt(sigma.var1*sigma.var2)
                 )
        return(out)
    }

    ## ** compute correlation and its iid
    e.cor <- FUN.cor(object.param)
    nabla <- numDeriv::jacobian(FUN.cor, object.param)
    colnames(nabla) <- name.param
    rownames(nabla) <- names(e.cor)
    
    if(robust){
        iid.cor <- object.iid %*% t(nabla)
    }else{
        iid.cor <- NULL
    }
    vcov.cor <- nabla %*% object.vcov %*% t(nabla)

    ## ** assemble
    out <- as.data.frame(matrix(NA, nrow = 6, ncol = 7, 
                                dimnames = list(c("variance 1","variance 2","direct covariance","total covariance","direct correlation","total correlation"),
                                                c("variable","estimate","se","lower","upper","null","p.value")                                  
                                                )))
    attr(out,"iid") <- iid.cor
    attr(out,"vcov") <- vcov.cor

    out$variable <- c(var1,var2,rep(paste0("(",var1,",",var2,")"),4))
    out$estimate <- e.cor
    out$se <- sqrt(diag(vcov.cor))
    out$lower <- out$estimate + qnorm(alpha/2) * out$se
    out$upper <- out$estimate + qnorm(1-alpha/2) * out$se
    out$null <- c(NA,NA,rep(null,4))
    out$p.value <- 2*(1-pnorm(abs(out$estimate-out$null)/out$se))

    ## ** export
    return(out)
}



##----------------------------------------------------------------------
### lvm2cor.R ends here
