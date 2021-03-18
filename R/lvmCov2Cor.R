### lvm2cor.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 17 2021 (20:28) 
## Version: 
## Last-Updated: mar 18 2021 (13:38) 
##           By: Brice Ozenne
##     Update #: 55
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
##' @param FUN the function used to compute the correlation. Alternative to the arguments var1 and var2.
##' @param FUN.args names of the coefficients used in FUN.
##'
##' @details WARNING by default, this function will "just" convert a covariance coefficient into a correlation coefficient,
##' i.e. divide it by the standard deviation of the corresponding variables.
##' It will **not** output the "global" correlation when the covariance is modeled through several parameters in the LVM (e.g. through a Latent Variable plus a specific covariance parameter).
##' 
##' @return A data.frame containing in the fourth row the estimated correlation coefficient, standard error, confidence intervals, and p.value.
##' @examples
##' #### Simulate some data ####
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
##' #### Fit the lvm ####
##' m1 <-lvm(c(PEQ_poslife,PEQ_posself,PEQ_posmood,PEQ_possoc,PEQ_posbehav)~lv.peq,
##'          c(MEQ_mystical,MEQ_mood) ~ 1*lv.meq,
##'          c(MEQ_timespace,MEQ_ineffability) ~ lv.meq)
##' covariance(m1) <- lv.peq ~ lv.meq
##' covariance(m1) <- MEQ_timespace~MEQ_ineffability
##' latent(m1) <- ~lv.peq + lv.meq
##' e <- estimate(m1, d)
##'
##' #### Delta method for the correlation between the LVs ####
##' ## approximated true value
##' c("var.meq" = var(GS$lv.meq), "var.peq" = var(GS$lv.peq),
##'   "cov" = cov(GS$lv.meq,GS$lv.peq), "cor" = cor(GS$lv.meq,GS$lv.peq))
##' 
##' ## using lava
##' estimate(e, function(x){
##' a <- x["lv.meq~~lv.meq"]
##' b <- x["lv.peq~~lv.peq"]
##' c <- x["lv.peq~~lv.meq"]
##'   c(var.meq = a,
##'    var.peq = b,
##'    cov = c,
##'    cor = c/sqrt(a*b))
##' })
##'
##' ## Delta method via lvmCov2Cor
##' lvmCov2Cor(e, var1 = "lv.meq", var2 = "lv.peq")
##' lvmCov2Cor(e, FUN = function(x){x["lv.peq~~lv.meq"]/sqrt(x["lv.peq~~lv.peq"]*x["lv.meq~~lv.meq"])})
##' 
##' #### Delta method for the correlation between endogenous ####
##' ## approximated true value
##' c("var.meq" = var(GS$MEQ_timespace), "var.peq" = var(GS$MEQ_ineffability),
##'   "cov" = cov(GS$MEQ_timespace,GS$MEQ_ineffability),
##'   "cor" = cor(GS$MEQ_timespace,GS$MEQ_ineffability))
##' 
##' ## using lava (partial correlation)
##' estimate(e, function(x){
##' a <- x["MEQ_timespace~~MEQ_timespace"]
##' b <- x["MEQ_ineffability~~MEQ_ineffability"]
##' c <- x["MEQ_timespace~~MEQ_ineffability"]
##'   c(var.meq = a,
##'    var.peq = b,
##'    cov = c,
##'    cor = c/sqrt(a*b))
##' })
##'
##' ## Delta method via lvmCov2Cor (partial correlation)
##' lvmCov2Cor(e, var1 = "MEQ_timespace", var2 = "MEQ_ineffability")
##' 
##' ## ERROR as no (direct) covariance parameter
##' ## lvmCov2Cor(e, var1 = "PEQ_poslife", var2 = "PEQ_posself")
##' 
##' ## Delta method via lvmCov2Cor (full correlation)
##' lvmCov2Cor(e, FUN = function(x){
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
##' return((sigma12+lambda1*lambda2*tau)/sqrt(Sigma1*Sigma2))
##' })
##' 


## * lvmCov2Cor (code)
##' @rdname lvm2cor
##' @export
lvmCov2Cor <- function(object, var1, var2, null = 0, level = 0.95, FUN = NULL, FUN.args = NULL){
    require(numDeriv)
    beta <- coef(object)
    phi <- iid(object)

    ## ** check arguments
    alpha <- 1-level
    if(is.null(FUN)){
    
        name.var1 <- paste0(var1,lava.options()$symbols[2],var1)
        name.var2 <- paste0(var2,lava.options()$symbols[2],var2)
        name.cov <- c(paste0(var1,lava.options()$symbols[2],var2),
                      paste0(var2,lava.options()$symbols[2],var1))
        name.cov <- name.cov[name.cov %in% names(beta)]
        test <- c(name.var1,name.var2,name.cov) %in% names(beta)
        if(name.var1 %in% names(beta) == FALSE){
            stop("Argument \'var1\' incorrect: no coefficient \"",var1,"\" could be find. \n")
        }
        if(name.var2 %in% names(beta) == FALSE){
            stop("Argument \'var2\' incorrect: no coefficient \"",var2,"\" could be find. \n")
        }
        if(name.var2 == name.var1){
            stop("Argument \'var2\' incorrect: should not be the same name as argument \'var1\'. \n")
        }
        if("cor" %in% c(name.var1,name.var2)){
            stop("No coefficient should be called \"cor\". This name is used internally. \n")
        }
        if(length(name.cov)==0){
            stop("Could not find a covariance parameter. \n")
        }

        FUN <- function(x){
            return(as.double(x[name.cov]/sqrt(prod(x[c(name.var1,name.var2)]))))
        }
        FUN.args <- c(name.var1,name.var2,name.cov)
    }else{
        if(!is.null(FUN.args)){
            if(any(FUN.args %in% beta == FALSE)){
                stop("Incorrect argument \'FUN.args\': should indicate coefficients from the argument \'object\'. \n",
                     "Coefficient(s) not found: \"",paste0(FUN.args[FUN.args %in% beta == FALSE], collapse ="\" \""),"\"\n")
            }
            if(any(duplicated(FUN.args))){
                stop("Incorrect argument \'FUN.args\': should not contain duplicated values. \n")
            }
        }
    }
        

    ## ** extract elements from object
    if(!is.null(FUN.args)){
        beta <- beta[FUN.args]
        phi <- phi[,FUN.args,drop=FALSE]
    }
    ## if(robust==FALSE){
    ##     vec.se <- sqrt(diag(crossprod(phi))/diag(vcov(object)))
    ##     k <- diag(1/vec.se)
    ##     phi <- phi %*% k
    ##     ## does not keep the correlations the same
    ##     ## crossprod(phi)-vcov(object)
    ## }
    ## ** compute correlation and its iid
    rho <- as.double(FUN(beta))
    if(is.na(rho)){
        stop("Could not compute the correlation coefficient using FUN. \n")
    }
    nabla <- numDeriv::jacobian(FUN, beta)
    phi.rho <- phi %*% t(nabla)
    colnames(phi.rho) <- "rho"
    if(is.null(FUN.args)){
        FUN.args <- names(beta)[abs(nabla)>1e-10 ]
        beta <- beta[FUN.args]
        phi <- phi[,FUN.args,drop=FALSE]
    }
    out <- as.data.frame(matrix(NA, nrow = length(FUN.args)+1, ncol=6, 
                                dimnames = list(c(FUN.args,"cor"),
                                                c("estimate","se","lower","upper","null","p.value")                                  
                                                )))
    attr(out,"iid") <- cbind(phi, phi.rho)

    ## ** assemble
    out$estimate <- c(beta, rho = rho)
    out$se <- sqrt(diag(crossprod(attr(out,"iid"))))
    out$lower <- out$estimate + qnorm(alpha/2) * out$se
    out$upper <- out$estimate + qnorm(1-alpha/2) * out$se
    out$null <- c(rep(NA,length(beta)),null)
    out$p.value <- 2*(1-pnorm(abs(out$estimate-out$null)/out$se))

    ## ** export
    return(out)
}



##----------------------------------------------------------------------
### lvm2cor.R ends here
