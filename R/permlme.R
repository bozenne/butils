### permlme.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 28 2021 (13:19) 
## Version: 
## Last-Updated: Jun  7 2021 (16:39) 
##           By: Brice Ozenne
##     Update #: 14
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * permlme (documentation)
##' @title Permutation test for fixed effects of a linear mixed model.
##' @description Permutation test for a single fixed effects of a linear mixed model.
##' This is essentially a re-implementation of the function \code{predictmeans::perlmer} limited to fixed effects and random intercept models.
##' Should be faster than the original function and output the test statistic relative to each permutation.
##' @name permlme
##' 
##' @param lme0 [lme] model under the null.
##' @param lme1 [lme] model under the alternative.
##' @param data [data.frame] dataset used to fit the models. Normally not need as the data is extracted from the models.
##' @param seed [integer] specify a random number generator seed, for reproducible results.
##' @param statistic [character] type of test statistic: Wald test and/or likelihood ratio test.
##' @param perm.X [logical] should the covariates be also permuted?
##' @param perm.index.perm [logical] should the index used to permute of the normalized residuals be output?
##' @param nperm [integer] number of permutations.
##' @param cpus [cpus] number of cpus used to perform in parallel the permutations.
##' @param trace [logical] should a progress bar displayed to follow the progress of the permutations?
##' 
##' @references Oliver E. Lee and Thomas M. Braun (2012), Permutation Tests for Random Effects in Linear Mixed Models. Biometrics, Journal 68(2).
##'
##' @return a list containing \itemize{
##' \item LRT: Likelihood ratio test for the coefficient. The p-value corresponding to the permutation test is denoted \code{p-value.perm}.
##' \item Wald: Wald test for the coefficient. The p-value corresponding to the permutation test is denoted \code{p-value.perm}.
##' \item n.perm: number of successful permutations
##' \item index.perm: index used to permute of the normalized residuals 
##' }

## * permlme (examples)
##' @rdname permlme
##' @examples
##' n.perm <- 1000
##' 
##' ## ** simulate data
##' library(lava)
##' library(data.table)
##' m <- lvm(Y1~1*eta+effect*conc,
##'          Y2~1*eta+effect*conc,
##'          Y3~1*eta+effect*conc,
##'          Y4~1*eta+effect*conc,
##'          Y5~1*eta+effect*conc)
##' latent(m) <- ~eta
##' transform(m, Id~eta) <- function(x,...){1:NROW(x)}
##'
##' set.seed(1)
##' dtW <- data.table::data.table(lava::sim(m, p = c("effect" = 0.1), n = 50, latent = FALSE))
##' dtL <- data.table::melt(dtW, id.vars = c("Id","conc"))
##' setkeyv(dtL,"Id")
##'
##' ## ** "homemade" permutation test
##' library(nlme)
##' e.lmeH0 <- nlme::lme(value ~ variable, random =~1|Id, data = dtL)
##' e.lmeH1 <- nlme::lme(value ~ variable + conc, random=~1|Id, data = dtL)
##' e.permlme <- permlme(e.lmeH0, e.lmeH1, perm.X = TRUE, nperm = n.perm, seed = 10, statistic = "LRT")
##' e.permlme ## p=0.2197802
##'
##' ## ** permutation test via the function permlmer
##' library(lme4)
##' library(predictmeans)
##' e.lmerH0 <- lme4::lmer(value ~ variable + (1|Id), data = dtL, REML = FALSE)
##' e.lmerH1 <- lme4::lmer(value ~ variable + conc + (1|Id), data = dtL, REML = FALSE)
##' e.permlmer <- predictmeans::permlmer(e.lmerH0, e.lmerH1, nperm = n.perm, ncore = 1, seed = 10)
##' e.permlmer
##'
##' ## ** compare values
##' e.permlme$LRT[2,"p.value.perm"]-e.permlmer[["Perm-p"]][2] ## same!!!
##' 
##' ## ** compare timing
##' library(microbenchmark)
##' microbenchmark(manual0 = permlme(e.lmeH0, e.lmeH1, nperm = 100, seed = 10, trace = FALSE),
##'                manual = permlme(e.lmeH0, e.lmeH1, nperm = 100, seed = 10, trace = FALSE, statistic = c("LRT","Wald")),
##'                package = predictmeans::permlmer(e.lmerH0, e.lmerH1, nperm = 100, ncore = 1, seed = 10),
##'                times = 10)
##' 
##' ## Unit: milliseconds (n=5)
##' ##    expr       min        lq      mean    median       uq       max neval cld
##' ## manual0  689.9769  695.0086  752.4042  749.6891  770.365  918.2436    10 a  
##' ##  manual 1963.7686 2053.5608 2111.0911 2135.7464 2147.678 2215.7608    10  b 
##' ## package 3433.9876 3456.7219 3591.3811 3641.0649 3662.663 3688.3372    10   c
##'
##' ## Unit: milliseconds
##' ## expr       min        lq     mean   median       uq      max neval cld
##' ## manual0  912.5592  939.0943 1117.806 1110.615 1278.572 1378.842    10 a
##' ## manual 2646.7307 3291.1398 3368.645 3352.098 3466.784 3829.754    10  b
##' ## package 4022.4287 4479.0607 4795.183 4530.408 4903.027 6810.617    10   c
##'
##' ## ** multiple testing
##'
##' ## *** unadjusted p-value
##' p.value <- e.permlme$LRT[["p.value.perm"]][2]
##'
##' ## *** perfect correlation between tests
##' e.permlme2 <- permlme(e.lmeH0, e.lmeH1, perm.X = TRUE, nperm = n.perm, seed = 10, statistic = "LRT")
##' e.permMax <- pmax(e.permlme$stat.perm[,"LRT"],e.permlme2$stat.perm[,"LRT"])
##' e.obsMax <- pmax(e.permlme$LRT[["L.Ratio"]][2],e.permlme2$LRT[["L.Ratio"]][2])
##' (sum(e.permMax > e.obsMax)+1)/(n.perm+1) - p.value ## should be 0
##' 
##' ## *** perfect independence between tests
##' set.seed(11)
##' dtW2 <- data.table::data.table(lava::sim(m, p = c("effect" = 0.1), n = 50, latent = FALSE))
##' dtL2 <- data.table::melt(dtW2, id.vars = c("Id","conc"))
##' setkeyv(dtL2,"Id")
##'
##' e.lmeH0.bis <- nlme::lme(value ~ variable, random =~1|Id, data = dtL2)
##' e.lmeH1.bis <- nlme::lme(value ~ variable + conc, random=~1|Id, data = dtL2)
##' e.permlme2 <- permlme(e.lmeH0.bis, e.lmeH1.bis, perm.X = TRUE, nperm = n.perm, seed = 10, statistic = "LRT")
##' p.value2 <- e.permlme2$LRT[["p.value.perm"]][2]
##' 
##' e.permMax <- pmax(e.permlme$stat.perm[,"LRT"],e.permlme2$stat.perm[,"LRT"])
##' e.obsMax <- pmax(e.permlme$LRT[["L.Ratio"]][2],e.permlme2$LRT[["L.Ratio"]][2])
##' (sum(e.permMax > e.obsMax)+1)/(n.perm+1) - 2*min(p.value ,p.value2) ## should be approximately true
##'
##' ## *** partially correlated tests
##' e.lmeH0.bis <- nlme::lme(value ~ variable, random =~1|Id, data = rbind(dtL,dtL2))
##' e.lmeH1.bis <- nlme::lme(value ~ variable + conc, random=~1|Id, data = rbind(dtL,dtL2))
##' e.permlme2 <- permlme(e.lmeH0.bis, e.lmeH1.bis, perm.X = TRUE, nperm = n.perm, seed = 10, statistic = "LRT")
##' p.value2 <- e.permlme2$LRT[["p.value.perm"]][2]
##'
##' e.permMax <- pmax(e.permlme$stat.perm[,"LRT"],e.permlme2$stat.perm[,"LRT"])
##' e.obsMax <- pmax(e.permlme$LRT[["L.Ratio"]][2],e.permlme2$LRT[["L.Ratio"]][2])
##' (sum(e.permMax > e.obsMax)+1)/(n.perm+1)/min(p.value ,p.value2) ## should be between 1 and 2
##' 
##' ## ** "homemade" permutation test in presence of missing values
##' set.seed(10)
##' dtL$valueMiss <- dtL$value
##' dtL$valueMiss[rbinom(NROW(dtL), size = 1, prob = 0.1)==1] <- NA
##' 
##' e.lmeH0 <- nlme::lme(valueMiss ~ variable, random =~1|Id, data = dtL, na.action = na.omit)
##' e.lmeH1 <- nlme::lme(valueMiss ~ variable + conc, random=~1|Id, data = dtL, na.action = na.omit)
##' e.permlme <- permlme(e.lmeH0, e.lmeH1, perm.X = TRUE, nperm = n.perm, seed = 10, statistic = "LRT")
##' e.permlme


## * permlme (code)
##' @rdname permlme 
##' @export
permlme <- function(lme0, lme1, data = NULL, seed = NULL, 
                    statistic = "Wald", perm.X = FALSE, return.index.perm = FALSE,
                    nperm = 1000, cpus = 1, trace = TRUE){ 

    statistic <- match.arg(statistic, c("LRT","Wald"), several.ok = TRUE)
    if(is.null(data)){data <- nlme::getData(lme1)}
    data <- as.data.frame(data)

    ## ** 0- normal test
    name.coef <- setdiff(names(fixef(lme1)),names(fixef(lme0)))
    if(length(name.coef)!=1 && "Wald" %in% statistic){
        stop("There should be exactly one parameter more in lme1 than in lme0. \n")
    }
    lme0.ML <- update(lme0, data = data, method = "ML")
    lme1.ML <- update(lme1, data = data, method = "ML")
    out <- list(call = match.call(),
                LRT = anova(lme0.ML, lme1.ML),
                Wald = summary(update(lme1, data = data, method = "REML"))$tTable[name.coef,,drop=FALSE],
                stat.perm = NULL,
                n.perm = c(LRT = NA, Wald = NA))
    
    ## ** 1- extract key quantities from input
    n.obs <- NROW(data) ## number of observations

    cluster <- getGroups(lme0) ## to which cluster (patient) each observation belongs
    U.cluster <- levels(cluster)
    n.cluster <- length(U.cluster)
    index.cluster <- lapply(U.cluster, function(iCluster){which(cluster==iCluster)}) ## used to restaure proper ordering after tapply

    Y <- getResponse(lme1) ## response
    name.Y <- all.vars(formula(lme1))[[1]] ## name of the response variable
    X0 <- model.matrix(formula(lme0),data) ## design matrix

    Omega0 <- getVarCov(lme0, type = "marginal", individuals = levels(cluster)) ## residual variance-covariance matrix
    beta0 <- fixef(lme0) ## regression coefficients

    ## ** 2- compute residuals
    Xbeta0 <- X0 %*% beta0
    residuals0 <- as.double(Y - X0 %*% beta0)

    ## ** 3- compute normalized residuals
    sqrtOmega0 <- lapply(Omega0,function(iOmega0){t(chol(iOmega0))})
    sqrtOmega0M1 <- lapply(sqrtOmega0,solve)

    residuals0N <- vector(length=n.obs, mode = "numeric")
    for(iCluster in 1:n.cluster){ ## iCluster <- 1
        residuals0N[index.cluster[[iCluster]]] <- sqrtOmega0M1[[iCluster]] %*% residuals0[index.cluster[[iCluster]]]
    }

    ## ** 4- estimate the distribution of the test statistics under the null
    warper <- function(iPerm){
        data.perm <- data

        ## permute residuals and fixed effects
        index.perm <- sample(1:n.obs)
        residuals0N.perm <- residuals0N[index.perm]

        ## rescale residuals
        for(iCluster in 1:n.cluster){ ## iCluster <- 1
            data.perm[[name.Y]][index.cluster[[iCluster]]] <- sqrtOmega0[[iCluster]] %*% residuals0N.perm[index.cluster[[iCluster]]]
        }
        ## add permuted fixed effects
        if(perm.X){
            data.perm[[name.Y]] <- data.perm[[name.Y]] + Xbeta0[index.perm,,drop=FALSE]
        }

        if("LRT" %in% statistic){
            lme0.permML <- try(update(lme0, data = data.perm, method = "ML"), silent = TRUE)
            lme1.permML <- try(update(lme1, data = data.perm, method = "ML"), silent = TRUE)
            if(inherits(lme0.permML,"try-error")||inherits(lme0.permML,"try-error")){
                LRT.stat <- NA
            }else{
                LRT.stat <- as.double(2*(logLik(lme1.permML)-logLik(lme0.permML)))
            }
        }else{
            LRT.stat <- NA
        }
        
        if("Wald" %in% statistic){
            lme1.permREML <- try(update(lme1, data = data.perm, method = "REML"), silent = TRUE)
            if(inherits(lme1.permREML,"try-error")){
                Wald.stat <- NA
            }else{
                Wald.stat <- summary(lme1.permREML)$tTable[name.coef,"t-value"]
            }
        }else{
            Wald.stat <- NA
        }
        iOut <- c(LRT = LRT.stat, Wald = Wald.stat)
        if(return.index.perm){
            attr(iOut,"index.perm") <- index.perm
        }
        return(iOut)
    }
    

    if(cpus==1 || is.null(cpus)){
        if(!is.null(seed)){set.seed(seed)}

        if(trace){
            require(pbapply)
            ls.perm <- pbapply::pblapply(1:nperm, function(i){warper(i)})
        }else{
            ls.perm <- lapply(1:nperm, function(i){warper(i)})
        }
    }else{
        cl <- snow::makeSOCKcluster(cpus)
        doSNOW::registerDoSNOW(cl)

        pb <- txtProgressBar(max = nperm, style=3)
        opts <- list(progress = function(n) setTxtProgressBar(pb, n))
        if(is.null(seed)){
            ls.perm <- foreach::`%dopar%`(
                                    foreach::foreach(i=1:nperm, .options.snow=opts, .packages = "nlme"), {
                                        warper(i)
                                    })
        }else{
            require(doRNG)
            set.seed(seed)
            ls.perm <- doRNG::`%dorng%`(
                                    foreach::foreach(i=1:nperm, .options.snow=opts, .packages = "nlme"), {
                                        warper(i)
                                    })
        }
        parallel::stopCluster(cl)
    }
    out$stat.perm <- do.call(rbind,ls.perm)
    if(return.index.perm){
        out$index.perm <- do.call(cbind,lapply(ls.perm,attr,"index.perm"))
    }
    
    ## ** 5- compare the observed statistic with its null distribution
    n.permLRT <- sum(!is.na(out$stat.perm[,"LRT"]))
    if(n.permLRT>0){
        out$LRT[["p.value.perm"]] <- c(NA,(sum(pmax(out$stat.perm[,"LRT"],0) >= out$LRT$L.Ratio[2], na.rm = TRUE) + 1)/(n.permLRT+1))
    }
    out$n.perm["LRT"]  <- n.permLRT

    n.permWald <- sum(!is.na(out$stat.perm[,"Wald"]))
    if(n.permWald>0){
        out$Wald <- data.frame(out$Wald, "p.value.perm" = (sum(abs(out$stat.perm[,"Wald"]) >= abs(out$Wald[,"t-value"]), na.rm = TRUE) + 1)/(n.permWald+1))
    }
        out$n.perm["Wald"]  <- n.permWald

    ## ** export
    class(out) <- append("permlme", class(out))
    return(out)
}

## * print.permlem
print.permlme <- function(x,...){
    if("p.value.perm" %in% names(x$LRT)){
        cat("Likelihood ratio test: \n")
        print(x$LRT)
    }
    if("p.value.perm" %in% names(x$Wald)){
        cat("Wald test: \n")
        print(x$Wald)
    }
    return(NULL)
}


##----------------------------------------------------------------------
### permlme.R ends here
