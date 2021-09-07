### matchPair.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 13 2019 (10:56) 
## Version: 
## Last-Updated: sep  7 2021 (11:35) 
##           By: Brice Ozenne
##     Update #: 85
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * matchPair (documentation)
##' @title McNemar Test for non-Independant Observations
##' @description Implementation of an extension of the McNemar test to non-independant observations proposed by Eliasziw et al. (1991).
##' @name matchPair
##' 
##' @param value [numeric vector] vector of binary values.
##' @param method [character vector] measurement method.
##' @param strata [character vector] index of the strata.
##' @param method.correlation [character] method used to compute the correlation, either \code{"full"} or \code{"discordant"}, from section 2 and 4 of Eliasziw et al. (1991), respectively.
##' Only relevant when \code{type} is \code{"correction"}.
##' @param type [character] approach used: \code{"correction"}, \code{"max-test"}, \code{"permutation"}.
##' @param n.perm [character] number of permutations used to compute the p-value. Only relevant when \code{type} is \code{"permutation"}.
##' @param id [character vector] index of the clusters. Only relevant when \code{type} is \code{"permutation"}.
##' @param ... additional arguments.
##'
##' @references Eliasziw M. and Donner A.. Application of the mcnemar test to non-independent matched pair data. Statistics in medicine, volume 10, 1981-1991 (1991)
##' 
##' @author Brice Ozenne

## * matchPair (examples)
##' @examples
##' n <- 100
##' 
##' set.seed(10)
##' n.obs <- rbinom(n, size = 2, prob = 1)+1
##' df <- data.frame(id = unlist(lapply(1:n, function(x){rep(x, n.obs[x])})),
##'            X = rbinom(sum(n.obs), size = 1, prob = 0.5),
##'            Y = rbinom(sum(n.obs), size = 1, prob = 0.5))
##' df$strata <- unlist(tapply(df$id, df$id, function(x){cumsum(duplicated(x))}))+1
##'
##' dfL <- data.table::melt(df, id.vars = c("id","strata"))
##' matchPair(value = dfL$value, method = dfL$variable, strata = dfL$strata, type = "correction")
##' matchPair(value = dfL$value, method = dfL$variable, strata = dfL$strata, type = "max-test")
##' matchPair(value = dfL$value, method = dfL$variable, strata = dfL$strata, id = df$id, type = "permutation")
##'
##' mcnemar.test(table(df$X,df$Y))

## * matchPair (code)
##' @rdname matchPair
##' @export
matchPair <- function(value, method, strata, type, ...){
    n <- length(value)
    Umethod <- unique(method)
    if(length(method)!=n){
        stop("length of \'method\' does not match length of \'value\' \n")
    }
    if(length(strata)!=n){
        stop("length of \'strata\' does not match length of \'value\' \n")
    }
    if(length(Umethod)!=2){
        stop("\'method\' must only take two possible values \n")
    }
    type <- match.arg(type, choices = c("correction","max-test","permutation"), several.ok = FALSE)

    out <- NULL
    if("correction" %in% type){
        out <- .matchPairCorr(value = value,
                              method = method,
                              strata = strata,
                              n = n,
                              Umethod = Umethod,
                              ...)
    }
    if("max-test" %in% type){
        out <- .matchPairMax(value = value,
                             method = method,
                             strata = strata,
                             n = n,
                             Umethod = Umethod,
                             ...)
    }
    if("permutation" %in% type){
        out <- .matchPairPerm(value = value,
                              method = method,
                              strata = strata,
                              Umethod = Umethod,
                              ...)
    }

    return(out)
}

## ** .matchPairCorr
##' @rdname matchPair
.matchPairCorr <- function(value, method, strata, method.correlation = "full", ...){

    dots <- list(...)
    n <- dots$n
    Umethod <- dots$Umethod
    
    ## ** compute sufficient statistics
    method.correlation <- match.arg(method.correlation, c("full","discordant"))
    table.strata <- tapply(1:n,strata, function(x){
        value.strata <- value[x]
        method.strata <- method[x]
        table(value.strata[method.strata==Umethod[1]],
              value.strata[method.strata==Umethod[2]])
    })

    K <- length(table.strata)
    a <- sapply(table.strata, function(x){x[1,1]})
    b <- sapply(table.strata, function(x){x[1,2]})
    c <- sapply(table.strata, function(x){x[2,1]})
    d <- sapply(table.strata, function(x){x[2,2]})

    n <- a+b+c+d
    S <- b+c ## number of discordant pairs
    Kd <- sum(S>=1) ## number of strata with discordant pairs 
    S.bar <- (1/Kd) * sum(S) ## average number of discordant pairs
    S0 <- S.bar - (sum(S-S.bar)^2 - (K-Kd)*S.bar^2)/(Kd*(Kd-1)*S.bar)

    ## ** estimate correlation
    if(method.correlation == "full"){
        
        ## section 4 in the article
        n.bar <- (1/K) * sum(n)
        n0 <- n.bar - sum((n-n.bar)^2)/(K*(K-1)*n.bar)
        P11 <- sum(a)/sum(n)
        P10 <- sum(b)/sum(n)
        P01 <- sum(c)/sum(n)
        P00 <- sum(d)/sum(n)

        BMS <- (1/K) * sum( ( (a - n * P11)^2 + (b - n * P10)^2 + (c - n*P01)^2 + (d - n*P00)^2 )/n )
        WMS <- (1/(K*(n.bar-1))) * sum( ( a*(n-a) + b*(n-b) + c*(n-c) + d*(n-d) )/n )

        rho.star <- (BMS - WMS)/(BMS + (n0-1)*WMS)
        rho <- 1 / (1 + P10*(1-rho.star)/rho.star + P01*(1-rho.star)/rho.star)
        
    }else if(method.correlation == "discordant"){
        ## section 2 in the article
        p.bar <- sum(b)/sum(S)

        BMS <- (1/Kd) * sum(((b-S*p.bar)^2/S)[S>=1])
        WMS <- (1/(Kd*(S.bar-1))) * sum((b*(S-b)/S)[S>=1])

        rho <- (BMS - WMS)/(BMS + (S0 - 1) * WMS)
        
    }

    ## ** original statistical test
    chi.mcNemar <- (sum(b)-sum(c))^2/(sum(b)+sum(c))
    chi.mcNemar.continuity <- (abs(sum(b)-sum(c))-1)^2/(sum(b)+sum(c))
    ## mcnemar.test(value[method==Umethod[1]],
    ##              value[method==Umethod[2]],
    ##              correct = FALSE)
    ## mcnemar.test(value[method==Umethod[1]],
    ##              value[method==Umethod[2]],
    ##              correct = TRUE)

    ## ** corrected statistical test
    nc <- S0 + Kd*(S.bar-S0)
    C.version1 <- 1 + rho * sum(S*(S-1))/sum(S) ## equals ## sum(S*(1+(S-1)*rho))/sum(S) ## 
    ## C.version2 <- 1 + (nc-1)*rho ## they are equal if the clusters have the same size.
    
    new.mcNemar <- chi.mcNemar/C.version1
    new.mcNemar.continuity <- chi.mcNemar.continuity/C.version1

    out <- data.frame(method = c("mcNemar","mcNemar.continuity","correct.mcNemar","corrected.mcNemar.continuity"),
                      statistic = c(chi.mcNemar,chi.mcNemar.continuity,new.mcNemar, new.mcNemar.continuity),
                      OR = sum(b)/sum(c),
                      correlation = c(0,0,rho,rho))
    out$df <- 1
    out$p.value <- 1-stats::pchisq(out$statistic, df = out$df)
    attr(out,"table") <- table.strata
    return(out)
}

## ** .matchPairMax
##' @rdname matchPair
.matchPairMax <- function(value, method, strata, ...){

    dots <- list(...)
    n <- dots$n
    Umethod <- dots$Umethod

    ls.strata <- tapply(1:n, strata, function(x){x} )

    ls.test <- lapply(ls.strata, function(x){
        iValue <- value[x]
        iMethod <- method[x]
        stats::mcnemar.test(iValue[iMethod == Umethod[1]],
                            iValue[iMethod == Umethod[2]])}
        )

    
    vec.value <- unlist(lapply(ls.test, "[[", "statistic"))
    vec.pvalue <- unlist(lapply(ls.test, "[[", "p.value"))
    vec.df <- 1
    vec.correctpvalue <- min(stats::p.adjust(vec.pvalue,method = "holm"))

    out <- data.frame(strata = unique(strata),
                      statistic = vec.value,
                      df = vec.df,
                      p.value = vec.pvalue,
                      adjusted.p.value = vec.correctpvalue)
    return(out)
}
## ** .matchPairPerm
##' @rdname matchPair
.matchPairPerm <- function(value, method, strata, id, n.perm = 1e3, ...){

    dots <- list(...)
    Umethod <- dots$Umethod

    dtL <- data.table::data.table(value = value,
                                  method = method,
                                  strata = strata,
                                  id = id)
    dtL[,c("index") := 1:.N]
    ls.indexX <- dtL[dtL$method == Umethod[1],list(list(.SD$index)),by = "id"][["V1"]]
    ls.indexY <- dtL[dtL$method == Umethod[2],list(list(.SD$index)),by = "id"][["V1"]]
    n.Id <- length(ls.indexX)

    ## e.test <- stats::mcnemar.test(dtL[dtL$method == Umethod[1],value],
    ## dtL[dtL$method == Umethod[2],value])
    ## ** statistic
    tableOrignal <- table(dtL[dtL$method == Umethod[1],.SD$value],
                          dtL[dtL$method == Umethod[2],.SD$value])
    statOrignal <- (abs(tableOrignal[1,2]-tableOrignal[2,1])-1)^2/(tableOrignal[1,2]+tableOrignal[2,1])
    OR <- tableOrignal[1,2]/tableOrignal[2,1]
    
    ## ** resampling    
    vec.statPerm <- sapply(1:n.perm, function(x){
        indicator.perm <- stats::rbinom(n.Id, prob = 0.5, size = 1)
            
        index.perm <- which(indicator.perm == 1)
        index.nperm <- which(indicator.perm == 0)

        new.indexX <- c(unlist(ls.indexX[index.nperm]),unlist(ls.indexY[index.perm]))
        new.indexY <- c(unlist(ls.indexY[index.nperm]),unlist(ls.indexX[index.perm]))

        tablePerm <- table(dtL$value[new.indexX],dtL$value[new.indexY])
        statPerm <- (abs(tablePerm[1,2]-tablePerm[2,1])-1)^2/(tablePerm[1,2]+tablePerm[2,1])
        return(statPerm)               
    })

    ## ** output
    out <- data.frame(method = "permutation", statistic = statOrignal, OR = OR, n.perm = n.perm,
                      p.value = mean(vec.statPerm > statOrignal) + 1/2 * mean(vec.statPerm == statOrignal))
    return(out)

}
######################################################################
### matchPair.R ends here
