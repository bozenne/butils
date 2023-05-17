### riskIDM.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 17 2023 (11:24) 
## Version: 
## Last-Updated: maj 17 2023 (17:05) 
##           By: Brice Ozenne
##     Update #: 23
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * riskIDM (example)
##' @title Estimate an Illness Death Model
##' @description Use multiple Cox models to fit an Illness Death Model from data in the wide format.
##' Output the fitted models and occupancy probabilities under different scenario (observed data, no transition to intermediate states, transition to a single intermediate state).
##'
##' @param formula [formula] A formula indicating baseline covariates on the right-hand side.
##' @param data [data.frame] dataset
##' @param time [numeric vector, >=0] times at which the occupancy probability should be estimated.
##' @param n.boot [interger, >=0] When strictly positive a non-parametric bootstrap is performed to quantify the uncertainty of the estimates.
##' The value then indicates the number of bootstrap samples.
##' @param var.id [character] name of the column containing the subject id, i.e. unique identifier for each line.
##' @param var.time [character vector of length 2] name of the columns containing the time variables, i.e. time at which each type of event happen (intermediate or absorbing).
##' If an intermediate event does not occur (e.g. no switch of treatment) then the time variable should be set to the end of follow-up time.
##' @param var.type [character vector of length 2] name of the columns containing event type indicator.
##' The first event type indicator can be categorical (multiple intermediate states) but the last one should be binary.
##' @param var.cov [character vecotr] optional baseline covariate values to be considered.
##' @param start.type [character] starting state. Deduced from \code{var.type} if left unspecified.
##' 
##' @rdname riskIDM
##' @examples
##'
##' library(survival)
##' library(mstate)
##' library(ggplot2)
##' library(riskRegression)
##'
##' #### data (remove ties) ###
##' data(ebmt3)
##' set.seed(10)
##' noise <- sort(rnorm(NROW(ebmt3$prtime),sd = 0.00001))
##' ebmt3$prtime <- ebmt3$prtime/365.25 + noise
##' ebmt3$rfstime <- ebmt3$rfstime/365.25 + noise
##'
##' #### PH without covariates ####
##' e.riskPH <- riskIDM(~1, data = ebmt3, PH = TRUE,
##'                     var.id = "id", n.boot = 100,
##'                     var.type = c("prstat", "rfsstat"),
##'                     var.time = c("prtime", "rfstime"))
##'
##' attr(e.riskPH, "model")
##' tail(e.riskPH[e.riskPH$scenario=="observed",])
##' plot(e.riskPH, type = "curve")
##' plot(e.riskPH, type = "curve", which = "all")
##' plot(e.riskPH, type = "stackplot")
##' plot(e.riskPH, type = "stackplot", which = "all")
##'
##' #### PH with covariates ####
##' eAdj.riskPH <- riskIDM(~age, data = ebmt3, PH = TRUE,
##'                        var.id = "id", 
##'                        var.type = c("prstat", "rfsstat"),
##'                        var.time = c("prtime", "rfstime"))
##' 
##' attr(eAdj.riskPH, "model")
##' tail(eAdj.riskPH[e.riskPH$scenario=="observed",])
##' plot(eAdj.riskPH, type = "curve")
##' plot(eAdj.riskPH, type = "stackplot")
##'
##' #### NPH without covariates ####
##' e.riskNPH <- riskIDM(~1, data = ebmt3, PH = FALSE,
##'                     var.id = "id", 
##'                     var.type = c("prstat", "rfsstat"),
##'                     var.time = c("prtime", "rfstime"))
##'
##' attr(e.riskNPH, "model")
##' tail(e.riskNPH[e.riskPH$scenario=="observed",])
##' plot(e.riskNPH, type = "curve")
##' plot(e.riskNPH, type = "stackplot")
##'
##' ## comparison with mstate
##' newdata.L <- data.frame(trans = c(1,2,3), strata = c(1,2,3))
##' tmat <- trans.illdeath(names = c("Tx", "PR", "RelDeath"))
##' 
##' msbmt <- msprep(time = c(NA, "prtime", "rfstime"),
##'                 status = c(NA, "prstat", "rfsstat"),
##'                 data = ebmt3, trans = tmat)
##'
##' ebmt.coxNPH <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data = msbmt)
##' ebmt.coxPH <- coxph(Surv(Tstart, Tstop, status) ~ from + strata(to), data = msbmt)
##' 
##' ebmt.msfitNPH <- msfit(ebmt.coxNPH, newdata = newdata.L, trans = tmat)
##' ebmt.probNPH <- probtrans(ebmt.msfitNPH, predt = 0)
##' plot(ebmt.probNPH, use.ggplot = TRUE)
##'
##' e.riskNPH.obs <- e.riskNPH[e.riskNPH$scenario=="observed",]
##' plot(ebmt.probNPH[[1]]$time, ebmt.probNPH[[1]]$pstate1, type = "l")
##' points(e.riskNPH.obs$time, e.riskNPH.obs$survival, type = "l", col = "red")
##' 
##' plot(ebmt.probNPH[[1]]$time, ebmt.probNPH[[1]]$pstate3, type = "l")
##' points(e.riskNPH.obs$time, e.riskNPH.obs$risk.2, type = "l", col = "red")
##' 
##' #### NPH with covariates ####
##' eAdj.riskNPH <- riskIDM(~age, data = ebmt3, PH = FALSE,
##'                        var.id = "id", 
##'                        var.type = c("prstat", "rfsstat"),
##'                        var.time = c("prtime", "rfstime"))
##' 
##' attr(eAdj.riskNPH, "model")
##' tail(eAdj.riskNPH[e.riskPH$scenario=="observed",])
##' plot(eAdj.riskNPH, type = "curve")
##' plot(eAdj.riskNPH, type = "stackplot")


## * riskIDM (code)
##' @name riskIDM
riskIDM <- function(formula, data, PH, time = NULL, n.boot = 0,
                    var.id, var.time, var.type, var.cov = NULL, start.type = NULL,
                    trace = TRUE){

    require(riskRegression)
    require(data.table)
    require(survival)
    require(prodlim)
    tol <- 1e-12
    original.time <- time

    ## ** normalize arguments
    ## formula
    if(is.list(formula)){
        if(length(formula)!=2){
            stop("When a list, argument \'formula\' should have length 2. \n")
        }
        if(is.null(names(formula))){
            names(formula) <- c("switch","outcome")
        }else if(any(sort(names(formula)) != sort(c("switch","outcome")))){
            stop("Names of argument \'formula\' should be \"switch\" or \"outcome\". \n")
        }
    }else if(inherits(formula,"formula")){
        formula <- list(switch = formula, outcome = formula)
    }

    ## var.cov
    if(is.null(var.cov)){
        var.cov <- union(all.vars(formula[[1]]),all.vars(formula[[2]]))
    }
    ## data
    data <- as.data.frame(data)
    rownames(data) <- NULL
    dataL <- reshapeIDM(data, var.id = var.id, var.time = var.time, var.type = var.type, var.cov = var.cov, start.type = start.type)
    n.obs <- NROW(data)

    ## states
    states <- attr(dataL,"states")
    n.switch <- length(states$switch)
    start.type <- states$censoring
    if(is.character(states$all)){
        states$name <- c(states$censoring, states$switch, states$outcome)        
        states$name.switch <- states$switch
    }else{
        states$name <- c("survival",paste0("risk.",states$switch), paste0("risk.",states$outcome))        
        states$name.switch <- paste0("state ", states$switch)
    }

    ## PH
    if(PH){
        cov.state <- "state.start"
    }else{
        cov.state <- "strata(state.start.num)"
    }

    ## ** prepare
    ls.formula <- vector(mode = "list", length = n.switch+1)
    for(iSwitch in 1:n.switch){ ## iSwitch <- 1
        ls.formula[[iSwitch]] <- as.formula(paste0("Surv(time.stop, state.stop.num==",iSwitch+1,") ~ ", strsplit(deparse(formula[[1]]), split = "~")[[1]][2]))
    }

    n.cov <- c(length(all.vars(formula[[1]])),length(all.vars(formula[[2]])))
    if(n.cov[2]==0){
        ls.formula[[n.switch+1]] <- as.formula(paste0("Surv(time.start, time.stop, state.stop.num==",n.switch+2,") ~ ",cov.state))
    }else{
        ls.formula[[n.switch+1]] <- as.formula(paste0("Surv(time.start, time.stop, state.stop.num==",n.switch+2,") ~ ",cov.state," + ",strsplit(deparse(formula[[2]]), split = "~")[[1]][2]))
    }
    ## always keep 0, the last observed time for each type of event, and all event times that are not censoring
    jump.time <- sort(unique(c(0,apply(data[,var.time],2,max),dataL$time.stop[dataL$state.stop!=0])))
    if(!is.null(time)){
        jump.timeR <- jump.time[jump.time<=max(time)]
    }else{
        jump.timeR <- jump.time
        time <- jump.time
    }
    if(length(jump.timeR)==0){
        stop("All requested timepoints are before the first event. \n")
    }
    scenarioNoSwitch <- paste0("no ",paste(states$name.switch, collapse = ", "))
    if(n.switch>1){
        scenarioSwitch <- sapply(states$name.switch, function(iS){paste0(iS," instead of ",paste(setdiff(states$name.switch,iS), collapse = ", "))})
    }else{
        scenarioSwitch <- NULL
    }
    scenarioAll <- c("observed",scenarioNoSwitch,scenarioSwitch)

    ## ** warper
    warper <- function(sample){
        ## *** data
        if(sample==0){
            iData <- data
            iDataL <- dataL
        }else{
            iData <- data[sample.int(NROW(data),replace = TRUE),,drop=FALSE]
            iData[[var.id]] <- 1:NROW(iData)
            iDataL <- reshapeIDM(iData, var.id = var.id, var.time = var.time, var.type = var.type, var.cov = var.cov, start.type = start.type)
        }
        iDataL0 <- iDataL[iDataL$state.start==states$censoring,]

        ## *** fit cox models
        e.switch <- lapply(1:n.switch, function(iSwitch){
            iFF.txt <- paste(deparse(ls.formula[[iSwitch]]), collapse = "") ## handle long formula that would be split on several lines
            e.outcome <- eval(parse(text=paste0("coxph(",iFF.txt,", data = iDataL0, y=TRUE, x=TRUE)")))
        })
        iFF.txt <- paste(deparse(ls.formula[[n.switch+1]]), collapse = "") ## handle long formula that would be split on several lines
        e.outcome <- eval(parse(text=paste0("coxph(",iFF.txt,", data = iDataL, y=TRUE, x=TRUE)")))

        ## *** extract hazards and evaluate risks
        if(n.cov[1]==0 && n.cov[2]==0){
            pred12 <- lapply(e.switch, function(iModel){predictCox(iModel, times = jump.timeR, type = "hazard")$hazard})
            lambda12 <- do.call(cbind,pred12)

            newdata.switch <- data.frame(state.start = factor(c(states$censoring, states$switch), levels = states$all),
                                         state.start.num = as.numeric(factor(c(states$censoring, states$switch), levels = states$all)))
            predX3 <- predictCox(e.outcome, times = jump.timeR, newdata = newdata.switch, type = "hazard")
            lambda13 <- predX3$hazard[1,]
            lambda23 <- t(predX3$hazard[-1,,drop=FALSE])
            out <- rbind(cbind(index.time = 1:length(jump.timeR),
                               .hazard2risk(jump.timeR, hazard12 = lambda12, hazard13 = lambda13, hazard23 = lambda23, states = states),
                               scenario = scenarioAll[1]),
                         cbind(index.time = 1:length(jump.timeR),
                               .hazard2risk(jump.timeR, hazard12 = 0*lambda12, hazard13 = lambda13, hazard23 = lambda23, states = states),
                               scenario = scenarioNoSwitch)
                         )
            if(n.switch>1){
                for(iS in 1:n.switch){
                    iLambda12 <- matrix(0, nrow = length(jump.timeR), ncol = n.switch)
                    iLambda12[,iS] <- rowSums(lambda12)
                    out <- rbind(out,
                                 cbind(index.time = 1:length(jump.timeR),
                                       .hazard2risk(jump.timeR, hazard12 = iLambda12, hazard13 = lambda13, hazard23 = lambda23, states = states),
                                       scenario = scenarioSwitch[iS])
                                 )
                }
            }

        }else{
            ## covariates categories
            iGrid.cov <- iData[,var.cov,drop=FALSE]
            iGridE.cov <- interaction(iGrid.cov, sep = "..X..")
            iTest.Ucov <- !duplicated(iGridE.cov)
            iId.Ucov <- iData[[var.id]][iTest.Ucov]
            iWeight.Ucov <- sapply(iGridE.cov[iTest.Ucov], function(iiUcov){mean(iiUcov==iGridE.cov)})
            iN.id.Ucov <- length(iId.Ucov)

            ## extract hazards
            iDataL0.red <- iDataL0[match(iId.Ucov, iDataL0[[var.id]]),,drop=FALSE]            
            pred12 <- do.call(rbind,
                              lapply(e.switch, function(iModel){predictCox(iModel, newdata = iDataL0.red, times = jump.timeR, type = "hazard")$hazard})
                              )
            ID.pred12 <- do.call(rbind,lapply(e.switch, function(iModel){iDataL0.red[var.id]}))
            indexID.pred12 <- by(1:NROW(ID.pred12),ID.pred12,identity)

            newdata.switch <- do.call(rbind,lapply(c(states$censoring, states$switch), function(iLevel){
                data.frame(state.start = factor(iLevel, levels = states$all),
                           state.start.num = as.numeric(factor(iLevel, levels = states$all)),
                           iDataL0.red[,c(var.id,var.cov),drop=FALSE])
            }))
            predX3 <- predictCox(e.outcome, times = jump.timeR, newdata = newdata.switch, type = "hazard")$hazard
            indexID.predX3 <- by(1:NROW(newdata.switch$id),newdata.switch$id,identity)

            ## prepare output
            iTemplate <- as.data.frame(matrix(0, nrow = length(jump.timeR), ncol = 2+1+n.switch+1,
                                              dimnames = list(NULL, c("index.time","time",states$name))))
            iTemplate$index.time <- 1:length(jump.timeR)
            iTemplate$time <- jump.timeR
            ls.out <- stats::setNames(lapply(scenarioAll, function(iSc){
                cbind(iTemplate, scenario = iSc)
            }), scenarioAll)

            ## hazard to risk
            for(iID in 1:iN.id.Ucov){ ## iID <- 1
                iLambda12 <- t(pred12[indexID.pred12[[iID]],,drop=FALSE])
                iLambda13 <- predX3[indexID.predX3[[iID]][1],]
                iLambda23 <- t(predX3[indexID.predX3[[iID]][-1],,drop=FALSE])
                iiWeight.Uvcov <- iWeight.Ucov[iID]

                ## evaluate risks
                iOut.observed <- .hazard2risk(jump.timeR, hazard12 = iLambda12, hazard13 = iLambda13, hazard23 = iLambda23, states = states)
                iOut.censoring <- .hazard2risk(jump.timeR, hazard12 = 0*iLambda12, hazard13 = iLambda13, hazard23 = iLambda23, states = states)

                ## update average risk
                ls.out[[scenarioAll[1]]][,states$name] <- ls.out[[scenarioAll[1]]][,states$name] + iOut.observed[,states$name] * iiWeight.Uvcov
                ls.out[[scenarioNoSwitch]][,states$name] <- ls.out[[scenarioNoSwitch]][,states$name] + iOut.censoring[,states$name] * iiWeight.Uvcov
                
                if(n.switch>1){
                    for(iS in 1:n.switch){
                        iSwitch <- states$switch[[iS]]

                        iiLambda12 <- matrix(0, nrow = length(jump.timeR), ncol = n.switch)
                        iiLambda12[,iS] <- rowSums(iLambda12)
                        
                        iOut.state <- .hazard2risk(jump.timeR, hazard12 = iiLambda12, hazard13 = iLambda13, hazard23 = iLambda23, states = states)
                        
                        ls.out[[scenarioSwitch[iS]]][,states$name] <- ls.out[[scenarioSwitch[iS]]][,states$name] + iOut.state[,states$name] * iWeight.Ucov
                    }
                }                         
            }
            out <- do.call(rbind,ls.out)
            rownames(out)<- NULL
        }

        ## *** export
        attr(out,"model") <- c(setNames(e.switch,states$switch),
                               setNames(list(e.outcome),states$outcome))
        return(out)
    }
        
    ## ** evaluate risks
    out <- warper(0)
    model <- attr(out,"model")

    if(n.boot>0){
        if(trace){
            require(pbapply)
            ls.out <- pbapply::pblapply(1:n.boot, warper)
        }else{
            ls.out <- lapply(1:n.boot, warper)
        }
        dt.out <- data.table::as.data.table(do.call(rbind,ls.out))
        se.out <- dt.out[,lapply(.SD, sd, na.rm = TRUE), by = c("time","index.time","scenario")]
        lower.out <- dt.out[,lapply(.SD, quantile, prob=0.025, na.rm = TRUE), by = c("time","index.time","scenario")]
        upper.out <- dt.out[,lapply(.SD, quantile, prob=0.975, na.rm = TRUE), by = c("time","index.time","scenario")]

        names(se.out)[match(states$name,names(se.out))] <- paste0(states$name,".se")
        names(lower.out)[match(states$name,names(lower.out))] <- paste0(states$name,".lower")
        names(upper.out)[match(states$name,names(upper.out))] <- paste0(states$name,".upper")

        out <- merge(out,se.out, by = c("index.time","time","scenario"))
        out <- merge(out,lower.out, by = c("index.time","time","scenario"))
        out <- merge(out,upper.out, by = c("index.time","time","scenario"))
    }

    if(!is.null(original.time)){
        out <- do.call(rbind,by(out,out$scenario, function(iOut){
            iiOut <- iOut[prodlim::sindex(jump.time = iOut[,"time"], eval.time = time),,drop=FALSE]
            iiOut$time <- time
            return(iiOut)
        }))
    }
    

    ## ** export
    attr(out,"model") <- model
    attr(out,"states") <- states
    attr(out,"jump.time") <- jump.time
    attr(out,"tol") <- tol
    attr(out,"n.boot") <- n.boot
    attr(out,"PH") <- PH
    class(out) <- append("riskIDM",class(out))
    return(out)
}

## * reshapeIDM (documentation)
##' @title Wide to Long Format For Illness Death Model
##' @description Convert a dataset from the wide to the long format when considering 1 starting state, 1 or many irreversible, exclusive, intermediate states, and 1 absorbing state.
##' @rdname reshapeIDM
##'
##' @param data [data.frame] data set in the wide format
##' @param var.id [character] name of the column containing the subject id, i.e. unique identifier for each line.
##' @param var.time [character vector of length 2] name of the columns containing the time variables, i.e. time at which each type of event happen (intermediate or absorbing).
##' If an intermediate event does not occur (e.g. no switch of treatment) then the time variable should be set to the end of follow-up time.
##' @param var.type [character vector of length 2] name of the columns containing event type indicator.
##' The first event type indicator can be categorical (multiple intermediate states) but the last one should be binary.
##' @param var.cov [character vecotr] optional baseline covariate values to be considered.
##' @param start.type [character] starting state. Deduced from \code{var.type} if left unspecified.
##'
##' @details Argument \code{var.time} and \code{var.type} should refer to the same states, first the intermediate states, and then the absorbing state.
##'
##' @return A data.frame with class \code{"dataIMD"}.
##'
##' @examples
##' #### generate data ####
##' set.seed(10)
##' n <- 1000 ## sample size (half)
##' tau <- 12 ## max follow-up time
##' Tevent <- rexp(2*n, rate = 1/10)
##' Tswitch <- c(runif(n, min = 1, max = 6), rep(Inf, n))
##' Cswitch <- sample.int(2, size = 2*n, replace = 2)
##' index.OC <- which(Tevent[1:n]>Tswitch[1:n])
##' Tevent[index.OC] <- Tswitch[index.OC] + rexp(length(index.OC), rate = c(1/5,1/2.5)[Cswitch[index.OC]])
##'
##' #### wide format ####
##' df.W <- data.frame(id = 1:(2*n),
##'                    gender = 0:1,
##'                    time.event = pmin(Tevent,tau),
##'                    time.switch = pmin(Tevent,Tswitch,tau),
##'                    switch = ifelse(Tswitch<pmin(Tevent,tau), Cswitch, 0),
##'                    event = as.numeric(Tevent <= tau))
##' head(df.W)
##' 
##' #### long format ####
##' df.L <- reshapeIDM(df.W,
##'                    var.id = "id",
##'                    var.type = c("switch","event"),
##'                    var.time = c("time.switch","time.event"),
##'                    var.cov = "gender")
##' head(df.L)
##'
##' #### in mstate ####
##' if(require(mstate)){
##' tmat <- transMat(x = list(c(2, 3, 4), c(4), c(4), c()), names = c("NoOC","IUC","OC", "Depression"))
##' tmat
##' }

## * reshapeIDM (code)
##' @name reshapeIDM
reshapeIDM <- function(data, var.id, var.time, var.type, var.cov = NULL, start.type = NULL){

    ## ** normalize arguments
    data <- as.data.frame(data)
    if(length(var.time)!=length(var.type)){
        stop("Arguments \'var.time\' and \'var.type\' should have the same length. \n")
    }
    if(length(var.type)!=2){
        stop("Arguments \'var.type\' should have length 2. \n")
    }
    if(var.id %in% names(data) == FALSE){
        stop("Argument \'var.id\' does not correspond to a column in argument \'data\' \n")
    }
    if(any(var.time %in% names(data) == FALSE)){
        stop("Argument \'var.time\' does not correspond to columns in argument \'data\' \n")
    }
    if(any(var.type %in% names(data) == FALSE)){
        stop("Argument \'var.type\' does not correspond to columns in argument \'data\' \n")
    }
    if(!is.null(var.cov) && any(var.cov %in% names(data) == FALSE)){
        stop("Argument \'var.cov\' does not correspond to columns in argument \'data\' \n")
    }
    if(is.numeric(data[[var.type[1]]]) && (any(data[[var.type[1]]] %% 1>0) || any(data[[var.type[1]]]<0))){
        stop("Argument \'var.type[1]\' should refer to a factor or a positive integer variable. \n")
    }else if(is.character(data[[var.type[1]]])){
        stop("Argument \'var.type[1]\' should refer to a factor or a positive integer variable. \n")
    }
    if(is.numeric(data[[var.type[2]]]) && any(data[[var.type[2]]] %in% 0:1 == FALSE)){
        stop("Argument \'var.type[2]\' should refer to a factor or a binary variable. \n")
    }else if(is.character(data[[var.type[2]]])){
        stop("Argument \'var.type[2]\' should refer to a factor or a binary variable. \n")
    }
    if(length(unique(data[[var.type[2]]]))%in% 1:2 == FALSE){
        stop("Argument \'var.type[2]\' should take one or two unique values. \n")
    }
    if(length(var.cov)==0){var.cov <- NULL}

    ## ** move to long format
    if(is.factor(data[[var.type[1]]])){
        if(is.null(start.type)){
            if(length(levels(data[[var.type[2]]]))==2){
                start.type <- levels(data[[var.type[2]]])[1]
            }else if(length(levels(data[[var.type[1]]]))>1){
                start.type <- levels(data[[var.type[1]]])[1]
            }
        }
        level.type1 <- setdiff(levels(data[[var.type[1]]]),start.type)
        level.type2 <- setdiff(levels(data[[var.type[2]]]),start.type)
        level.all <- c(start.type, level.type1, level.type2)

    }else{
        if(!is.null(start.type)){
            if(start.type!=0){
                stop("Argument \'start.type\' must be 0 when using a numeric variable to encode types. \n")
            }
        }else{
            start.type <- 0
        }
        level.type1 <- setdiff(unique(data[[var.type[1]]]),0)
        level.type2 <- length(level.type1)+1
        level.all <- c(start.type, level.type1, level.type2)
        data[[var.type[2]]][data[[var.type[2]]]==1] <- level.type2
        
    }

    data.0 <- data[data[[var.type[1]]]==start.type,]
    data.1 <- data[data[[var.type[1]]]!=start.type,]


    out.1 <- data.frame(setNames(list(data.0[[var.id]]),var.id),
                        time.start = 0,
                        time.stop = data.0[[var.time[2]]],
                        state.start = start.type,
                        state.stop = data.0[[var.type[2]]])
    out.2 <- data.frame(setNames(list(data.1[[var.id]]),var.id),
                        time.start = 0,
                        time.stop = data.1[[var.time[1]]],
                        state.start = start.type,
                        state.stop = data.1[[var.type[1]]])
    out.3 <- data.frame(setNames(list(data.1[[var.id]]),var.id),
                        time.start = data.1[[var.time[1]]],
                        time.stop = data.1[[var.time[2]]],
                        state.start = data.1[[var.type[1]]],
                        state.stop = data.1[[var.type[2]]])
    if(!is.null(var.cov)){
        out.1 <- cbind(out.1,data.0[var.cov])
        out.2 <- cbind(out.2,data.1[var.cov])
        out.3 <- cbind(out.3,data.1[var.cov])
    }

    out <- rbind(out.1,out.2,out.3)
    out <- out[order(out[[var.id]],out$time.start),]
    rownames(out) <- NULL
    attr(out, "var") <- list(id = var.id,
                             cov = var.cov)
    attr(out, "states") <- list(all = level.all,
                                censoring = setdiff(level.all, c(level.type1,level.type2)),
                                switch = level.type1,
                                outcome = level.type2)

    ## ** export
    out$state.start.num <- as.numeric(factor(out$state.start, levels = level.all))
    out$state.start <- as.factor(out$state.start)
    out$state.stop.num <- as.numeric(factor(out$state.stop, levels = level.all))
    class(out) <- append("dataIDM",class(out))
    return(out)
}

## * plot.riskIDM
plot.riskIDM <- function(x, ...){
    require(ggplot2)
    out <- autoplot.riskIDM(x, ...)
    print(out)
    return(invisible(out))
}

## * autoplot.riskIDM (documentation)
##' @title Graphical display for For Illness Death Model
##' @description Diplay state occupation probability of an Illness Death Model
##' @rdname autoplot.riskIDM
##'
##' @param object [riskIDM] output of the \code{riskIDM} function.
##' @param type [character] type of plot:
##' occupation probability of a given state under various scenario (\code{"curve"})
##' or state occupation probability stacked for a specific scenario (\code{"stackplot"}).
##' @param which the scenario (\code{type="stackplot"}) or state (\code{type="curve"}) for which the occupation probabilities are displayed.
##' Can also be \code{"all"} to display all possibilites using facets.
##' @param ci [logical] should pointwise confidence intervals be displayed. Only available for \code{type="curve"}
##' when a non-parametric bootstrap has been performed when running  the \code{riskIDM} function.
##' @param linewidth [numeric, >0] thickness of the line used to display the occupation probabilities.
##' Only relevant for \code{type="curve"}.
##' @param ci.alpha [numeric, 0-1] transparency parameter for the pointwise confidence intervals.
##' @param breaks [numeric vector, 0-1] labels used for the y-axis
##' @param ... not used
##' 
##' @return A ggplot2 object.
##'
##' @examples
##' library(survival)
##' library(mstate)
##' library(ggplot2)
##' library(riskRegression)
##'
##' #### data (remove ties) ###
##' data(ebmt3)
##' set.seed(10)
##' noise <- sort(rnorm(NROW(ebmt3$prtime),sd = 0.00001))
##' ebmt3$prtime <- ebmt3$prtime/365.25 + noise
##' ebmt3$rfstime <- ebmt3$rfstime/365.25 + noise
##'
##' #### fit IDM ####
##' e.riskPH <- riskIDM(~1, data = ebmt3, PH = FALSE,
##'                     var.id = "id", 
##'                     var.type = c("prstat", "rfsstat"),
##'                     var.time = c("prtime", "rfstime"))
##'
##' e.riskPH2 <- riskIDM(~1, data = ebmt3, PH = FALSE, n.boot = 100,
##'                     var.id = "id",  
##'                     var.type = c("prstat", "rfsstat"),
##'                     var.time = c("prtime", "rfstime"))
##' 
##' #### graphical display ####
##' plot(e.riskPH, type = "stackplot")
##' plot(e.riskPH, type = "stackplot", which = "all")
##' 
##' plot(e.riskPH, type = "curve")
##' plot(e.riskPH2, type = "curve", ci = TRUE)
##' plot(e.riskPH2, type = "curve", ci = TRUE, which = "all")


## * autoplot.riskIDM (code)
##' @name autoplot.riskIDM
autoplot.riskIDM <- function(object, type = "stackplot", which = NULL, ci = TRUE,
                             linewidth = 2, ci.alpha = 0.2, breaks = seq(0,1,by=0.1), ...){

    ## ** normalize arguments
    type <- match.arg(type, c("stackplot","curve"))

    dots <- list(...)
    if (length(dots) > 0) {
        stop("Unknown argument(s) '", paste(names(dots), collapse = "' '"), "'. \n")
    }

    ## ** graphical display
    states <- attr(object, "states")
    jump.time <- attr(object, "jump.time")
    tol <- attr(object, "tol")
    n.boot <- attr(object, "n.boot")
    if(n.boot==0){
        ci <- FALSE
    }
    
    if(type == "curve"){
        ## prepare
        if(is.null(which)){
            state <- utils::tail(states$name,1)
            label.state <- utils::tail(states$all,1)
            label.y <- paste0("Occupancy probability for state \"",label.state,"\"") 
        }else{
            if(identical(which,"all")){
                which <- states$all
            }
            label.state <- match.arg(as.character(which), as.character(states$all), several.ok = TRUE)
            state <- states$name[label.state == as.character(states$all)]
            if(is.numeric(states$all)){
                label.state2 <- paste0("state ", label.state)
            }else{
                label.state2 <- label.state
            }
            if(length(state)==1){
                label.y <- paste0("Occupancy probability for state \"",label.state,"\"") 
            }else{
                label.y <- "Occupancy probability"
            }
        }
        ## reshape data to long format
        outL <- reshape2::melt(object[,c("index.time","time","scenario",state)],
                               id.vars = c("index.time","time","scenario"),
                               value.name = "estimate",
                               variable.name = "state")
        outL$state <- factor(outL$state, levels = state, labels = label.state2)
        if(ci){
            outL.ci <- cbind(reshape2::melt(object[,c("index.time","time","scenario",paste0(state,".lower"))],
                                            id.vars = c("index.time","time","scenario"),
                                            value.name = "lower", variable.name = "state"),
                             upper = reshape2::melt(object[,c("index.time","time","scenario",paste0(state,".upper"))],
                                                    id.vars = c("index.time","time","scenario"),
                                                    value.name = "upper", variable.name = "state")$upper)
            outL.ci$state <- factor(outL.ci$state, levels = paste0(state,".lower"), labels = label.state2)
        }
        ## reshape2::dcast(object[object$time>11.9,c("time","survival","scenario")],
        ##                 formula = scenario ~ time, value.var = "survival")

        ## graphical display
        gg <- ggplot2::ggplot(outL, ggplot2::aes(x = time, group = scenario))
        if(n.boot>0 & ci){
            gg <- gg + ggplot2::geom_ribbon(data = outL.ci, ggplot2::aes(ymin = lower, ymax = upper, fill = scenario), alpha = ci.alpha)
        }
        gg <- gg + ggplot2::geom_line(linewidth = linewidth, ggplot2::aes(y = estimate, color = scenario))
        gg <- gg + ggplot2::scale_y_continuous(breaks = breaks, labels=scales::percent)
        gg <- gg + ggplot2::labs(fill = "Scenario", color = "Scenario", y = label.y)
        if(length(label.state)>1){
            gg <- gg + ggplot2::facet_wrap(~state)
        }

    }else if(type == "stackplot"){
        valid.scenario <- unique(object$scenario)
        if(is.null(which)){
            which <- "observed"
        }else if(identical(which, "all")){
            which <- valid.scenario
        }
        scenario <- match.arg(which, valid.scenario, several.ok = TRUE)
        if(length(scenario)==1){
            label.y <- paste0("Occupancy probability for scenario \"",scenario,"\"") 
        }else{
            label.y <- "Occupancy probability"
        }
        out1 <- object[object$scenario %in% scenario,c("index.time","time","scenario",states$name)]
        out1L <- reshape2::melt(out1, id.vars = c("index.time","time","scenario"))
        if(is.character(states$all)){
            out1L$variable <- factor(out1L$variable, levels = rev(c(states$censoring, states$switch, states$outcome)))
        }else{
            out1L$variable <- factor(out1L$variable,
                                     levels = rev(states$name),
                                     labels = rev(states$all))
        }
        out1L.after <- as.data.frame(out1L)
        out1L.after$time <- c(jump.time[-1]-tol,tail(jump.time,1)-tol)[out1L.after$index.time]

        gg <- ggplot2::ggplot(rbind(out1L,out1L.after), ggplot2::aes(x = time, y = value, group = variable, fill = variable))
        gg <- gg + ggplot2::geom_area() + ggplot2::labs(fill = "State", y = label.y)
        gg <- gg + ggplot2::scale_y_continuous(breaks = breaks, labels=scales::percent)
        gg <- gg + ggplot2::coord_cartesian(ylim = c(0,1))
        if(length(scenario)>1){
            gg <- gg + ggplot2::facet_wrap(~scenario)
        }
    }

    ## ** export
    return(gg)
}

## * .hazard2risk
.hazard2risk <- function(time, hazard12, hazard13, hazard23, states, prodlim = TRUE){

    n.time <- length(time)
    n.switch <- NCOL(hazard12)

    cumhazard12 <- colCumSum(hazard12)
    cumhazard13 <- cumsum(hazard13)
    cumhazard23 <- colCumSum(hazard23)

    ## ** pre-compute
    if(prodlim){
        S11 <- cumprod(1-rowSums(hazard12)-hazard13)
    }else{
        S11 <- exp(- rowSums(cumhazard12) - cumhazard13)
    }
    S11m <- c(1,S11[-n.time])
    S11m.hazard12 <- riskRegression::colMultiply_cpp(hazard12, S11m)

    if(prodlim){
        ## ## fast implementation
        C1.cumhazard23 <- apply(1-hazard23,2,cumprod)
        S11m.hazard12.C1.cumhazard23.inv <- riskRegression::colCumSum(S11m.hazard12/C1.cumhazard23)
        S12 <- S11m.hazard12.C1.cumhazard23.inv * C1.cumhazard23

        ## ## ## slow but explicit implementation
        ## range(S12 - do.call(rbind,lapply(1:n.time, function(iTau){ ## iTau <- 10
        ##     if(iTau==1){return(0)} ## both hazard cannot be simulataneously non-0 at time 1
        ##     iScale <- matrix(apply(1-hazard23[1:iTau,,drop=FALSE],2,prod), nrow = iTau, ncol = n.switch, byrow = TRUE)
        ##     iS22 <- iScale/apply(1-hazard23[1:iTau,,drop=FALSE],2,cumprod)            
        ##     iFactor <- matrix(S11m[1:iTau], nrow = iTau, ncol = n.switch, byrow = FALSE)
        ##     return(colSums(iFactor * hazard12[1:iTau,,drop=FALSE] * iS22[1:iTau,,drop=FALSE]))
        ## })))

    }else{ ## old slow version still necessary when prodlim == FALSE
        ## ## fast implementation
        e.cumhazard23 <- exp(cumhazard23)
        S11m.hazard12.e.cumhazard23 <-  riskRegression::colCumSum(S11m.hazard12 * e.cumhazard23)
        S12 <- S11m.hazard12.e.cumhazard23 / e.cumhazard23

        ## ## slow but explicit implementation
        ## range(S12 - do.call(rbind,lapply(1:n.time, function(iTau){ ## iTau <- 10
        ##     if(iTau==1){return(0)} ## both hazard cannot be simulataneously non-0 at time 1
        ##     iCenter <- matrix(cumhazard23[iTau,,drop=FALSE], nrow = iTau, ncol = n.switch, byrow = TRUE)
        ##     iS22 <- exp(-iCenter+cumhazard23[1:iTau,,drop=FALSE])
        ##     iFactor <- matrix(S11m[1:iTau], nrow = iTau, ncol = n.switch, byrow = FALSE)
        ##     return(colSums(iFactor * hazard12[1:iTau,,drop=FALSE] * iS22[1:iTau,,drop=FALSE]))
        ## })))
    }
    S12m <- rbind(rep(0,n.switch),S12[-n.time,,drop=FALSE])
    
    ## ** scenario
    out <- data.frame(time,
                      S11,
                      S12,
                      cumsum(hazard13 * S11m) + rowSums(colCumSum(hazard23 * S12m))
                      )
    
    ## ** export
    names(out) <- c("time",states$name)
    return(out)
}

## * .commonString
## find the common consecutive part between two strings
## copied from https://stackoverflow.com/questions/35381180/identify-a-common-pattern
## .commonString("aaachang2","aaabbb")
## .commonString("aaa235change2","aaachangebbb")
## .commonString("abcdef","xyz")
.commonString <- function(string1, string2){

    A <- strsplit(string1, "")[[1]]
    B <- strsplit(string2, "")[[1]]

    L <- matrix(0, length(A), length(B))
    ones <- which(outer(A, B, "=="), arr.ind = TRUE)
    ones <- ones[order(ones[, 1]), ,drop=FALSE]
    if(NROW(ones)==0){return(NA)}
    for(i in 1:NROW(ones)) {
        v <- ones[i, , drop = FALSE]
        L[v] <- ifelse(any(v == 1), 1, L[v - 1] + 1)
    }
    out <- paste0(A[(-max(L) + 1):0 + which(L == max(L), arr.ind = TRUE)[1]], collapse = "")
    return(out)
}


##----------------------------------------------------------------------
### riskIDM.R ends here
