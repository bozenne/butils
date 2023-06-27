### riskIDM.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj 17 2023 (11:24) 
## Version: 
## Last-Updated: jun 27 2023 (16:05) 
##           By: Brice Ozenne
##     Update #: 521
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * riskIDM (documentation)
##' @title Estimate an Illness Death Model
##' @description Use multiple Cox models to fit an Illness Death Model from data in the wide format.
##' Output the fitted models and occupancy probabilities under different scenario (observed data, no transition to intermediate states, transition to a single intermediate state).
##' @rdname riskIDM
##'
##' @param formula [formula] A formula indicating baseline covariates on the right-hand side.
##' @param data [data.frame] dataset
##' @param PH [logical] should the ratio between the hazard with and without exposure be time independent?
##' @param time [numeric vector, >=0] times at which the occupancy probability should be estimated.
##' @param intervention [list of matrix] list where each element is a matrix used to deduce the intervention hazards by premultiplying the estimated hazard.
##' @param n.boot [interger, >=0] When strictly positive a non-parametric bootstrap is performed to quantify the uncertainty of the estimates.
##' The value then indicates the number of bootstrap samples.
##' @param level [numeric, 0.1] Confidence level for the confidence interval.
##' @param cl [interger, >=0] A cluster object created by \code{makeCluster}, or an integer to indicate number of child-processes (integer values are ignored on Windows) for parallel evaluations
##' Passed to \code{pblapply}. Ignored when \code{trace=FALSE}.
##' @param var.id [character] name of the column containing the subject id, i.e. unique identifier for each line.
##' @param var.time [character vector of length 2] name of the columns containing the time variables, i.e. time at which each type of event happen (intermediate or absorbing).
##' If an intermediate event does not occur (e.g. no switch of treatment) then the time variable should be set to the end of follow-up time.
##' @param var.type [character vector of length 2] name of the columns containing event type indicator.
##' The first event type indicator can be categorical (multiple intermediate states) but the last one should be binary.
##' @param start.type [character] starting state. Deduced from \code{var.type} if left unspecified.
##' @param keep.indiv [logical] should covariate specific occupancy probabilities be output?
##' @param trace [logical] should a progress bar be used to display the execution of the resampling procedure?

## * riskIDM (examples)
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
##' summary(e.riskPH)
##' model.tables(e.riskPH)
##' confint(e.riskPH, time = 1:10)
##' confint(e.riskPH, time = 2:11)
##' coef(e.riskPH, time = 1:10)
##' model.tables(e.riskPH, contrast = "all", time = c(1,2))
##' confint(e.riskPH, contrast = "all", time = c(1,2))
##' coef(e.riskPH, contrast = "all", time = c(1,2))
##' plot(e.riskPH, by = "scenario")
##' plot(e.riskPH, by = "scenario", scenario = "all")
##' plot(e.riskPH, by = "state")
##' plot(e.riskPH, by = "state", state = "all")
##' plot(e.riskPH, by = "contrast")
##'
##' #### PH with covariates ####
##' dt.ebmt3 <- as.data.table(ebmt3)
##' eAdj.riskPH <- riskIDM(~age, data = dt.ebmt3, PH = TRUE,
##'                        var.id = "id", 
##'                        var.type = c("prstat", "rfsstat"),
##'                        var.time = c("prtime", "rfstime"))
##' 
##' coef(eAdj.riskPH, time = 1:10)
##' plot(eAdj.riskPH, by = "scenario")
##' plot(eAdj.riskPH, by = "state")
##'
##' #### NPH without covariates ####
##' e.riskNPH <- riskIDM(~1, data = ebmt3, PH = FALSE,
##'                     var.id = "id", 
##'                     var.type = c("prstat", "rfsstat"),
##'                     var.time = c("prtime", "rfstime"))
##'
##' plot(e.riskNPH)
##' plot(e.riskNPH, by = "state")
##' plot(e.riskNPH, by = "contrast")
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
##' e.survNPH.obs <- confint(e.riskNPH, state = "survival")[scenario == "observed",]
##' plot(ebmt.probNPH[[1]]$time, ebmt.probNPH[[1]]$pstate1, type = "l")
##' points(e.survNPH.obs$time, e.survNPH.obs$estimate, type = "l", col = "red")
##' 
##' e.riskNPH.obs <- confint(e.riskNPH, state = "risk.2")[scenario == "observed",]
##' plot(ebmt.probNPH[[1]]$time, ebmt.probNPH[[1]]$pstate3, type = "l")
##' points(e.riskNPH.obs$time, e.riskNPH.obs$estimate, type = "l", col = "red")
##' 
##' #### NPH with covariates ####
##' eAdj.riskNPH <- riskIDM(~age, data = ebmt3, PH = FALSE,
##'                        var.id = "id", 
##'                        var.type = c("prstat", "rfsstat"),
##'                        var.time = c("prtime", "rfstime"))
##' 
##' summary(eAdj.riskNPH)
##' coef(eAdj.riskNPH)
##' plot(eAdj.riskNPH, by = "scenario")
##' plot(eAdj.riskNPH, by = "state", state = "all")
##'
##' #### multiple exposures ####
##' set.seed(10)
##' n <- 1000 ## sample size (half)
##' tau <- 12 ## max follow-up time
##' Tevent <- rexp(2*n, rate = 1/10)
##' Tswitch <- c(runif(n, min = 1, max = 6), rep(Inf, n))
##' Cswitch <- sample.int(2, size = 2*n, replace = 2)
##' index.OC <- which(Tevent[1:n]>Tswitch[1:n])
##' Tevent[index.OC] <- Tswitch[index.OC] + rexp(length(index.OC), rate = c(1/5,1/2.5)[Cswitch[index.OC]])
##' 
##' df.W <- data.frame(id = 1:(2*n),
##'                    gender = 0:1,
##'                    time.event = pmin(Tevent,tau),
##'                    time.switch = pmin(Tevent,Tswitch,tau),
##'                    switch = ifelse(Tswitch<pmin(Tevent,tau), Cswitch, 0), 
##'                    event = as.numeric(Tevent <= tau))
##' df.W$switch <- factor(df.W$switch, levels = 0:2, c(0,"OC","IUD"))
##' # df.W$switch <- factor(df.W$switch, levels = 0:2, c("H","OC","IUD"))
##' # df.W$event <- factor(df.W$event, levels = 0:1, c("H","MDD"))
##' eME.riskPH <- riskIDM(~1, data = df.W, PH = FALSE,
##'                     var.id = "id",
##'                     var.type = c("switch","event"),
##'                     var.time = c("time.switch","time.event"))
##'
##' summary(eME.riskPH)
##' confint(eME.riskPH, contrast = "all")
##' plot(eME.riskPH, by = "scenario")
##' plot(eME.riskPH, by = "scenario", scenario = "all")
##' plot(eME.riskPH, by = "state")
##' plot(eME.riskPH, by = "state", state = "all")
##' plot(eME.riskPH, by = "contrast", scenario = "all")
##' plot(eME.riskPH, by = "contrast", scenario = c("observed","no OC, IUD"))


## * riskIDM (code)
##' @rdname riskIDM
##' @export
riskIDM <- function(formula, data, PH, time = NULL, intervention = NULL,
                    var.id, var.time, var.type, start.type = NULL, 
                    n.boot = 0, level = 0.95, cl = NULL,
                    keep.indiv = FALSE, trace = TRUE){

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
    var.cov <- union(all.vars(formula[[1]]),all.vars(formula[[2]]))
    if(!is.null(var.cov)){
        if(keep.indiv && any(var.cov=="weight")){
            stop("Inconstency between argument \'data\' and argument \'keep.indiv\'. \n",
                 "\"weight\" is being used internally and should not be a column name. \n")
        }
        if(keep.indiv && any(var.cov=="signature")){
            stop("Inconstency between argument \'data\' and argument \'keep.indiv\'. \n",
                 "\"signature\" is being used internally and should not be a column name.")
        }
        if(is.character(keep.indiv)){
            sep.cov <- keep.indiv
        }else if(keep.indiv){
            sep.cov  <- "."
        }else{
            sep.cov <- NULL
        }
    }else{
        sep.cov <- NULL
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

    ## intervention
    name.lambda01 <- paste0("lambda.0",1:n.switch) ## healthy -> illness
    name.lambda02 <- paste0("lambda.0",n.switch+1) ## healthy -> death
    name.lambda12 <- paste0("lambda.",1:n.switch,n.switch+1) ## illness -> death
    name.lambda <- c(name.lambda01,name.lambda02,name.lambda12)
    n.lambda <- length(name.lambda)
        
    if(is.null(intervention)){
        scenarioNoSwitch <- paste0("no ",paste(states$name.switch, collapse = ", "))
        if(n.switch>1){
            scenarioSwitch <- sapply(states$name.switch, function(iS){paste0(iS," instead of ",paste(setdiff(states$name.switch,iS), collapse = ", "))})
        }else{
            scenarioSwitch <- NULL
        }
        scenarioAll <- unname(c("observed",scenarioNoSwitch,scenarioSwitch))
        n.intervention <- length(scenarioAll)

        intervention <- stats::setNames(lapply(1:n.intervention, function(iSc){
            matrix(0, nrow = n.lambda, ncol = n.lambda,
                   dimnames = list(name.lambda, name.lambda)
                   )
        }), scenarioAll)

        diag(intervention[["observed"]]) <- 1 ## no change (identity matrix)
        diag(intervention[[scenarioNoSwitch]]) <- c(rep(0,n.switch),rep(1,1+n.switch)) ## set transition to exposure (i.e. illness states) to 0
        if(!is.null(scenarioSwitch)){ ## transfer transition to other exposures to the one exposure (i.e. illness states 2,3,4 ---> illness state 1)
            for(iSc in 1:length(scenarioSwitch)){ ## iSc <- scenarioSwitch[1]
                diag(intervention[[scenarioSwitch[iSc]]]) <- c(rep(0,n.switch),rep(1,1+n.switch))
                intervention[[scenarioSwitch[iSc]]][iSc,] <- c(rep(1,n.switch),rep(0,1+n.switch))
            }
        }        
    }else{
        if(is.matrix(intervention)){
            intervention <- list("manual" = intervention)
        }
        if(!is.list(intervention)){
            stop("Argument \'intervention' should be a list. \n")
        }
        if(is.null(names(intervention))){
            stop("Argument \'intervention\' should be a named list. \n")
        }
        if(any(duplicated(names(intervention)))){
            stop("Argument \'intervention\' should be a named list with non-duplicated names. \n")
        }
        scenarioAll <- names(intervention)
        if(any(sapply(intervention,is.matrix)==FALSE)){
            stop("Argument \'intervention\' should be a list of matrices. \n")
        }
        if(any(sapply(intervention,dim)!=n.lambda)){
            stop("Argument \'intervention\' should be a list of matrices of size ",n.lambda," by ",n.lambda,". \n")
        }
        n.intervention <- length(scenarioAll)

        for(iI in 1:n.intervention){
            if(is.null(colnames(intervention[[iI]]))){
                colnames(intervention[[iI]]) <- name.lambda
            }else if(any(colnames(intervention[[iI]])!=name.lambda)){
                stop("Argument \'intervention\' should be a list of matrices with column names \"",paste0(name.lambda, collapse = "\", \""),"\". \n")
            }
            if(is.null(rownames(intervention[[iI]]))){
                rownames(intervention[[iI]]) <- name.lambda
            }else if(any(colnames(intervention[[iI]])!=name.lambda)){
                stop("Argument \'intervention\' should be a list of matrices with row names \"",paste0(name.lambda, collapse = "\", \""),"\". \n")
            }
        }
    }
    
    ## ** prepare
    ls.formula <- vector(mode = "list", length = n.switch+1)
    for(iSwitch in 1:n.switch){ ## iSwitch <- 1
        ls.formula[[iSwitch]] <- stats::update(formula[[1]], paste0("Surv(time.stop, state.stop.num==",iSwitch+1,") ~ ."))
    }

    n.cov <- c(length(all.vars(formula[[1]])),length(all.vars(formula[[2]])))
    if(n.cov[2]==0){
        ls.formula[[n.switch+1]] <- as.formula(paste0("Surv(time.start, time.stop, state.stop.num==",n.switch+2,") ~ ",cov.state))
    }else{
        ls.formula[[n.switch+1]] <- stats::update(formula[[2]], paste0("Surv(time.start, time.stop, state.stop.num==",n.switch+2,") ~ ",cov.state," + ."))
    }
    ## always keep 0, the last observed time for each type of event, and all event times that are not censoring
    
    jump.time <- unique(sort(c(0, ## first timepoint
                               tapply(dataL$time.stop,dataL$state.stop,max), ## last timepoint
                               dataL[as.character(dataL$state.start)!=as.character(dataL$state.stop),"time.stop"]))) ## time point for each change of state

    if(!is.null(time)){
        jump.timeR <- jump.time[jump.time<=max(time)]
    }else{
        jump.timeR <- jump.time
        time <- jump.time
    }
    if(length(jump.timeR)==0){
        stop("All requested timepoints are before the first event. \n")
    }
    
    ## ** warper
    warper <- function(sample, sep.cov){
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
        out <- NULL
        if(n.cov[1]==0 && n.cov[2]==0){
            pred01 <- lapply(e.switch, function(iModel){predictCox(iModel, times = jump.timeR, type = "hazard")$hazard})
            lambda01 <- do.call(cbind,pred01)

            newdata.switch <- data.frame(state.start = factor(c(states$censoring, states$switch), levels = states$all),
                                         state.start.num = as.numeric(factor(c(states$censoring, states$switch), levels = states$all)))
            predX3 <- predictCox(e.outcome, times = jump.timeR, newdata = newdata.switch, type = "hazard")
            lambda02 <- predX3$hazard[1,]
            lambda12 <- t(predX3$hazard[-1,,drop=FALSE])

            M.lambda <- cbind(lambda01,lambda02,lambda12)
            
            for(iSc in scenarioAll){ ## iSc <- scenarioAll[1]
                iLambda01 <- M.lambda %*% base::t(intervention[[iSc]][name.lambda01,,drop=FALSE])
                iLambda02 <- M.lambda %*% base::t(intervention[[iSc]][name.lambda02,,drop=FALSE])
                iLambda12 <- M.lambda %*% base::t(intervention[[iSc]][name.lambda12,,drop=FALSE])

                out <- rbind(out,
                             cbind(index.time = 1:length(jump.timeR),
                                   .hazard2risk(jump.timeR, hazard01 = iLambda01, hazard02 = iLambda02, hazard12 = iLambda12, states = states),
                                   scenario = iSc)
                             )
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
            pred01 <- do.call(rbind,
                              lapply(e.switch, function(iModel){predictCox(iModel, newdata = iDataL0.red, times = jump.timeR, type = "hazard")$hazard})
                              )
            iDataL0.red <- as.data.frame(iDataL0.red) ## ensures that predictCox did not transform it into a data.table
            ID.pred01 <- do.call(rbind,lapply(e.switch, function(iModel){iDataL0.red[var.id]}))
            indexID.pred01 <- by(1:NROW(ID.pred01),ID.pred01,identity)

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
            if(!is.null(sep.cov)){
                attr(ls.out,"indiv") <- cbind(iDataL0.red[,c(var.id,var.cov),drop=FALSE],
                                              weight = iWeight.Ucov,
                                              signature = as.character(interaction(iDataL0.red[,var.cov], sep = sep.cov)))
            }
                

            ## hazard to risk
            for(iID in 1:iN.id.Ucov){ ## iID <- 1
                iLambda01 <- t(pred01[indexID.pred01[[iID]],,drop=FALSE])
                iLambda02 <- predX3[indexID.predX3[[iID]][1],]
                iLambda12 <- t(predX3[indexID.predX3[[iID]][-1],,drop=FALSE])
                iM.lambda <- cbind(iLambda01, iLambda02, iLambda12)
                iiWeight.Uvcov <- iWeight.Ucov[iID]

                for(iSc in scenarioAll){ ## iSc <- scenarioAll[1]
                    iiLambda01 <- iM.lambda %*% base::t(intervention[[iSc]][name.lambda01,,drop=FALSE])
                    iiLambda02 <- iM.lambda %*% base::t(intervention[[iSc]][name.lambda02,,drop=FALSE])
                    iiLambda12 <- iM.lambda %*% base::t(intervention[[iSc]][name.lambda12,,drop=FALSE])

                    iOut <- .hazard2risk(jump.timeR, hazard01 = iiLambda01, hazard02 = iiLambda02, hazard12 = iiLambda12, states = states)
                    ls.out[[iSc]][,states$name] <- ls.out[[iSc]][,states$name] + iOut[,states$name] * iiWeight.Uvcov
                    if(!is.null(sep.cov)){
                        iSignature <- attr(ls.out,"indiv")$signature[iID]
                        attr(ls.out,iSignature) <- rbind(attr(ls.out,iSignature),
                                                         cbind(iOut, scenario = iSc, signature = iSignature)
                                                         )
                    }
                }
            }
            out <- do.call(rbind,ls.out)
            rownames(out)<- NULL

            if(!is.null(sep.cov)){
                attr(out,"indiv") <- do.call(rbind,lapply(attr(ls.out,"indiv")$signature, function(iSignature){attr(ls.out,iSignature)}))
                rownames(attr(out,"indiv")) <- NULL
            }
            
        }

        ## *** export
        attr(out,"model") <- c(setNames(e.switch,states$switch),
                               setNames(list(e.outcome),states$outcome))
        return(out)
    }
        
    ## ** evaluate risks
    ## point estimate
    res <- warper(0, sep.cov = sep.cov)
    ## store 
    out <- list(args = list(PH = PH, time = time, intervention = intervention, n.boot = n.boot, level = level,
                            var.id = var.id, var.time = var.time, var.type = var.type, start.type = start.type, indiv = keep.indiv),
                call = match.call(),
                data = dataL,
                jump.estimate = as.data.table(res[,c("index.time","time","scenario",states$name)]),
                jump.indiv = as.data.table(attr(res,"indiv")[,c("time","scenario","signature",states$name)]),
                jump.time = jump.time,
                model = attr(res,"model"),
                scenario = scenarioAll,
                states = states,
                tol = tol
                )
    ## move to long fomat
    out$jump.estimate <- data.table::melt(out$jump.estimate,
                                          id.vars = c("index.time","time","scenario"),
                                          measure.vars = states$name,
                                          value.name = "estimate",
                                          variable.name = "state")
    if(keep.indiv){
        out$jump.indiv <- data.table::melt(out$jump.indiv,
                                           id.vars = c("time","scenario","signature"),
                                           measure.vars = states$name,
                                           value.name = "estimate",
                                           variable.name = "state")
    }
    ## ** evaluate uncertainty
    if(n.boot>0){
        alpha <- 1-level
        if(trace){
            require(pbapply)
            ls.out <- pbapply::pblapply(1:n.boot, function(iBoot){
                iRes <- try(warper(iBoot, sep.cov = NULL))
                if(inherits(iBoot,"try-error")){
                    return(NULL)
                }else{
                    return(cbind(boot = iBoot, iRes))
                }
            }, cl = cl)
        }else{
            ls.out <- lapply(1:n.boot, function(iBoot){
                iRes <- try(warper(iBoot, sep.cov = NULL))
                if(inherits(iBoot,"try-error")){
                    return(NULL)
                }else{
                    return(cbind(boot = iBoot, iRes))
                }
            } = NULL)
        }
        out$boot <- data.table::as.data.table(do.call(rbind,ls.out)[,c("boot","index.time","time","scenario",states$name)])

        dt.bootmedian <- out$boot[,lapply(.SD, median, na.rm = TRUE), by = c("time","scenario"), .SDcols = states$name]
        dt.bootlower <- out$boot[,lapply(.SD, quantile, 1-alpha/2, na.rm = TRUE), by = c("time","scenario"), .SDcols = states$name]
        dt.bootupper <- out$boot[,lapply(.SD, quantile, alpha/2, na.rm = TRUE), by = c("time","scenario"), .SDcols = states$name]

        dt.boot <- data.table::melt(dt.bootmedian,
                                    id.vars = c("time","scenario"),
                                    measure.vars = states$name,
                                    value.name = "median",
                                    variable.name = "state")
        dt.boot$lower <- data.table::melt(dt.bootlower,
                                          id.vars = c("time","scenario"),
                                          measure.vars = states$name,
                                          value.name = "lower",
                                          variable.name = "state")$lower
        dt.boot$upper <- data.table::melt(dt.bootupper,
                                          id.vars = c("time","scenario"),
                                          measure.vars = states$name,
                                          value.name = "upper",
                                          variable.name = "state")$upper
        out$jump.estimate <- merge(x = out$jump.estimate, y = dt.boot, by = c("time","scenario","state"))
    }

    ## ** export
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
##' tau <- 01 ## max follow-up time
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
##' @rdname reshapeIDM
##' @export
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
    if(is.factor(data[[var.type[1]]]) && is.factor(data[[var.type[2]]])){
        if(levels(data[[var.type[1]]])[1]!=levels(data[[var.type[2]]])[1]){
            stop("Variables defined by \'var.type\' should have the same reference level. \n")
        }
        if(levels(data[[var.type[2]]])[2] %in% levels(data[[var.type[1]]])){
            stop("The second level of the second variable of \'var.type\' not match any level of the first variable of \'var.type\'. \n")
        }
    }
    if(!is.null(start.type) && is.factor(data[[var.type[1]]]) && start.type %in% levels(data[[var.type[1]]]) == FALSE){
        stop("Argument \'start.type\' does not match any level of the first variable of \'var.type\'. \n")
    }
    if(length(var.cov)==0){var.cov <- NULL}

    ## ** move to long format
    if(is.factor(data[[var.type[1]]]) && is.factor(data[[var.type[2]]])){
        if(is.null(start.type)){
            if(length(levels(data[[var.type[2]]]))==2){
                start.type <- levels(data[[var.type[2]]])[1]
            }else{
                start.type <- levels(data[[var.type[1]]])[1]
            }
        }
        level.type1 <- setdiff(levels(data[[var.type[1]]]),start.type)
        level.type2 <- setdiff(levels(data[[var.type[2]]]),start.type)
        level.all <- c(start.type, level.type1, level.type2)

    }else if(is.factor(data[[var.type[1]]])){
        if(is.null(start.type)){
            start.type <- levels(data[[var.type[1]]])[1]
        }
        level.type1 <- setdiff(levels(data[[var.type[1]]]),start.type)
        level.type2 <- length(level.type1)+1
        level.all <- c(start.type, level.type1, level.type2)

        data[[var.type[2]]] <- factor(data[[var.type[2]]], levels = 0:1, labels = c(start.type,length(level.type1)+1))
    }else if(is.factor(data[[var.type[2]]])){
        if(is.null(start.type)){
            if(length(levels(data[[var.type[2]]]))==2){
                start.type <- levels(data[[var.type[2]]])[1]
            }else{
                start.type <- 0
            }
        }
        level.type1 <- 1:(length(unique(data[[var.type[1]]]))-1)
        level.type2 <- setdiff(levels(data[[var.type[2]]]),start.type)
        level.all <- c(start.type, level.type1, level.type2)

        data[[var.type[1]]] <- factor(data[[var.type[1]]], levels = 0:length(level.type1), labels = c(start.type,level.type1))
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
    index.noswitch <- which(data[[var.time[1]]]==data[[var.time[2]]])
    if(any(data[index.noswitch,var.type[1]]!=start.type)){
        stop("Some patients had both the terminal event and switched to an intermediate state with a single time to event (",var.time[1],"=",var.time[2],"). \n",
             "Maybe something went wrong when identifying the reference state (here ",start.type,"). \n")
    }


    data.0 <- data[data[[var.type[1]]]==start.type,]
    data.1 <- data[data[[var.type[1]]]!=start.type,]
    if(NROW(data.0)>0){
        out.1 <- data.frame(setNames(list(data.0[[var.id]]),var.id),
                            time.start = 0,
                            time.stop = data.0[[var.time[2]]],
                            state.start = start.type,
                            state.stop = data.0[[var.type[2]]])
    }else{
        out.1 <- NULL
    }
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
        if(!is.null(out.1)){
            out.1 <- cbind(out.1,data.0[var.cov])
        }
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
    out$state.start <- factor(out$state.start, levels = level.all[-length(level.all)])
    out$state.stop.num <- as.numeric(factor(out$state.stop, levels = level.all))
    out$state.stop <- factor(out$state.stop, levels = level.all)
    class(out) <- append("dataIDM",class(out))
    
    return(out)
}

## * print.riskIDM
print.riskIDM <- function(x, ...){
    summary(x, short = TRUE)

    return(invisible(NULL))
}


summary.riskIDM <- function(object, short = FALSE, time = NULL, digits = c(2,3)){

    n.switch <- length(object$state$switch)
    if(n.switch==1){
        cat("     Illness Death Model with 1 intermediate state. \n")
    }else{
        cat("     Illness Death Model with ",n.switch," intermediate states \n\n",sep="")
    }
    cat(" - ",NROW(object$data)," observations from ",length(unique(object$data[[object$args$var.id]]))," individuals\n", sep = "")
    stoptable <- table(object$data[["state.stop"]])
    cat("   number of stops per state: ",paste(paste0(names(stoptable)," = ",stoptable), collapse =", "),"\n", sep = "")
    
    cat(" - ",length(object$jump.time)," timepoints: 0 to ",max(object$jump.time),"\n", sep = "")
    cat(" - ",length(object$scenario)," scenarios: \"",paste(object$scenario, collapse = "\"\n                \""),"\"\n", sep = "")
    if(length(unlist(lapply(object$model,coef)))==0){
        cat(" - non-parametric hazard estimators \n", sep = "")
    }else{
        cat(" - semi-parametric hazard estimator \n", sep = "")
        if(!short){
            
            ls.model <- stats::setNames(lapply(object$model,function(iModel){
                summary(iModel)$coefficient
            }), object$state$name[-1])
            ls.model <- ls.model[!sapply(ls.model,is.null)]
            print(ls.model)
        }
    }
    if(short){
        cat(" - estimated state occupancy (observed scenario)\n", sep = "")
        print(object$jump.estimate[scenario=="observed",.(min = 100*min(estimate, na.rm = TRUE),
                                                          median = 100*median(estimate, na.rm = TRUE),
                                                          max = 100*max(estimate, na.rm = TRUE)),by="state"],
              row.names = FALSE)
    }else{
        cat(" - estimated state occupancy under each scenario\n", sep = "")
        if(length(digits)==1){
            digits <- rep(digits,2)
        }
        if(is.null(time)){
            time <- c(quantile(object$jump.time, probs = c(0,0.25,0.5,0.75)),
                      object$jump.estimate[,.(NNA=sum(!is.na(estimate))), by = "time"][NNA>0,max(time)])
        }
        dt.state <- do.call(rbind,lapply(object$states$name, function(iState){
            model.tables(object, time = time, state = iState)
        }))
        dt.state$time <- round(dt.state$time, digits = digits[2])
        dt.state$index.time <- NULL
        dt.state$estimate <- round(100*dt.state$estimate, digits = digits[1])
        if(object$args$n.boot>0){
            dt.state$estimate <- paste0(dt.state$estimate," [",round(100*dt.state$lower, digits = digits[1]),";",round(100*dt.state$upper, digits = digits[1]),"]")
        }
        dtW.state <- data.table::dcast(data = dt.state, value.var = "estimate", formula = scenario+time ~ state)
        dtW.state[duplicated(dtW.state$scenario), scenario := ""]
        print(dtW.state, row.names = FALSE)
    }
    return(invisible(NULL))
}



## * model.tables.riskIDM (documentation)
##' @title Extract Probabilities From IDM
##' @description Extract occupancy probabilities, difference, or ratio between probabilities from an illness death model.
##' @rdname model.tables.riskIDM
##' 
##' @param object [riskIDM] output of the \code{riskIDM} function.
##' @param time [numeric vector] time at which the probabilities should be extracted. Can be \code{"unique"} to extract at jump times (i.e. non-duplicated risk values).
##' @param indiv [logical] should covariate specific probabilities be extracted?
##' @param state [character] state relative to which the occupancy probabilities should be extracted.
##' @param contrast [character vector] optional argument indicating scenario that are to be compared.
##' @param metric [character] how to compare the probabilities between two scenarios: \code{"difference"} (alternative - reference) or \code{"ratio"} (alternative / reference).
##' 

## * model.tables.riskIDM (code)
##' @rdname model.tables.riskIDM
##' @export
model.tables.riskIDM <- function(x, time = "unique", indiv = FALSE, state = utils::tail(x$states$name,1),
                                 contrast = NULL, metric = "difference", ...){

    ## ** normalize arguments
    dots <- list(...)
    if(length(dots) > 0) {
        stop("Unknown argument(s) '", paste(names(dots), collapse = "' '"), "'. \n")
    }

    x.args <- x$args
    if(!identical(indiv,FALSE)){
        if(x.args$indiv==FALSE){
            stop("Incorrect value for argument \'indiv\': individual occupancy probabilities have not been stored. \n",
                 "Consider setting the argument \'keep.indiv\' to TRUE when calling riskIDM. \n")
        }
        table <- data.table::copy(x$jump.indiv)
    }else{
        table <- data.table::copy(x$jump.estimate)
    }
    x.state <- x$state
    if(is.numeric(state)){
        state <- x.state$name[state]
    }
    state <- match.arg(state, x.state$name)

    name.scenario <- x$scenario
    if(!is.null(contrast)){
        if(!identical(indiv,FALSE)){
            stop("Cannot contrast individual occupancy probabilities. \n")
        }
        if(identical(contrast,"all")){
            contrast <- name.scenario
        }
        if(length(contrast)<2){
            stop("Argument \'contrast\' should contain at least two elements, e.g. ",name.scenario[1]," and ",name.scenario[2],". \n")
        }
        if(length(contrast)>length(name.scenario)){
            stop("Argument \'contrast\' should contain at most ",length(name.scenario)," elements. \n")
        }
        if(any(duplicated(contrast))){
            stop("Argument \'contrast\' should not contain duplicated values. \n")
        }
        contrast <- match.arg(contrast, name.scenario, several.ok = TRUE)
    }
    metric <- match.arg(metric, c("difference","ratio"))
        
    ## ** subset
    ## subset in two steps to avoid confusion between the argument (states) and the column names (since data.table is used)
    if(is.character(indiv)){
        index.keep <- intersect(which(table$state %in% state),which(table$signature %in% indiv))
    }else{
        index.keep <- which(table$state %in% state)
    }
    table <- table[index.keep]

    ## ** extract estimates
    Utime <- sort(unique(table$time))
    if(identical(time,"unique")){
        if(indiv){
            time <- unname(sort(unique(do.call(c,by(table, interaction(table$scenario,table$state,table$signature), function(iDF){
                iDF$time[!duplicated(iDF$estimate)]
            })))))
        }else{
            time <- unname(sort(unique(do.call(c,by(table, interaction(table$scenario,table$state), function(iDF){
                iDF$time[!duplicated(iDF$estimate)]
            })))))
        }
    }
    if(!is.null(time)){
        UindexTime <- prodlim::sindex(jump.time = Utime, eval.time = time)
        Utime.original <- Utime[UindexTime]
        if(indiv){
            out <- do.call(rbind,by(table, interaction(table$scenario,table$state,table$signature), function(iDF){
                iiDF <- iDF[match(Utime.original,iDF$time)]
                iiDF$time <- time
                return(iiDF)
            }))
        }else{
            out <- do.call(rbind,by(table, interaction(table$scenario,table$state), function(iDF){
                iiDF <- iDF[match(Utime.original,iDF$time)]
                iiDF$time <- time
                return(iiDF)
            }))
        }
    }else{
        out <- table
        time <- Utime
    }
    data.table::setcolorder(out, neworder = c("state","scenario","time","index.time"))
    out$scenario <- factor(out$scenario, levels = x$scenario)
    data.table::setkeyv(out, cols = c("scenario","time"))

    ## contrast between exposures
    if(!is.null(contrast)){

        allContrast <- combn(contrast,2)
        n.contrast <- NCOL(allContrast)

        out2 <- lapply(1:n.contrast, function(iC){ ## iC <- 1
            iContrast1 <- allContrast[1,iC]
            iContrast2 <- allContrast[2,iC]

            ## estimate
            if(metric == "difference"){
                iOut <- data.table::data.table(state =  state, time = time, alternative = iContrast2, reference = iContrast1,
                                               estimate = out[out$scenario == iContrast2,.SD$estimate] - out[out$scenario == iContrast1,.SD$estimate]
                                               )
            }else if(metric == "ratio"){
                iOut <- data.table::data.table(state =  state, time = time, alternative = iContrast2, reference = iContrast1,
                                               estimate = out[out$scenario == iContrast2,.SD$estimate] / out[out$scenario == iContrast1,.SD$estimate]
                                               )
                iOut$ratio[out[out$scenario == iContrast2,.SD$estimate]==0] <- 0
            }
            ## uncertainty
            if(x.args$n.boot>0){
                alpha <- 1-x.args$level
                iBoot1 <- x$boot[scenario == iContrast1,.SD,.SDcols = c("boot","time",state)]
                iBoot2 <- x$boot[scenario == iContrast2,.SD,.SDcols = c("boot","time",state)]
                iIndex.keep1 <- which(iBoot1$time %in% Utime.original)
                iIndex.keep2 <- which(iBoot2$time %in% Utime.original)
                
                if(metric == "difference"){
                    iBoot <- data.table::data.table(time = time,
                                        estimate = iBoot2[[state]][iIndex.keep2] - iBoot1[[state]][iIndex.keep1])
                }else if(metric == "ratio"){
                    iBoot <- data.table::data.table(time = time,
                                        estimate = ifelse(iBoot2[[state]][iIndex.keep2]==0,0,iBoot2[[state]][iIndex.keep2] / iBoot1[[state]][iIndex.keep1])
                                        )
                }
                iOut <- merge(x = iOut,
                              y = iBoot[,.(median = median(.SD$estimate, na.rm = TRUE),
                                           lower = quantile(.SD$estimate, probs = alpha/2, na.rm = TRUE),
                                           upper = quantile(.SD$estimate, probs = 1-alpha/2, na.rm = TRUE)),
                                        by="time"],
                              by = "time")
            }
            return(iOut)
        })
        out <- do.call(rbind,out2)
        data.table::setcolorder(out, neworder = c("state","alternative","reference","time"))
    }

    ## ** export
    return(out)
}

## * confint.riskIDM (code)
##' @export
confint.riskIDM <- function(object, ...){

    out <- model.tables(object, ...)
    out$state <- NULL
    if("index.time"  %in% names(out)){
        out$index.time <- NULL
    }
    if("median" %in% names(out)){
        out$median <- NULL
    }
    return(out)
}

## * coef.riskIDM (code)
##' @export
coef.riskIDM <- function(object, contrast = NULL, metric = "difference", ...){

    outAll <- model.tables(object, contrast = contrast, metric = metric, ...)

    Utime <- unique(outAll$time)
    n.time <- length(Utime)
    if(!is.null(contrast)){
        if(metric == "difference"){
            outAll$scenario <- paste0(outAll$alternative,"-",outAll$reference)
        }else if(metric == "ratio"){
            outAll$scenario <- paste0(outAll$alternative,"/",outAll$reference)
        }
    }
    Uscenario <- unique(outAll$scenario)
    n.scenario <- length(Uscenario)
    out <- matrix(outAll$estimate, nrow = n.time, ncol = n.scenario, byrow = FALSE,
                  dimnames = list(Utime, Uscenario))
    return(out)
}

## * model.frame.riskIDM
##' @export
model.frame.riskIDM <- function(formula, ...){
    formula$model
}

## * plot.riskIDM
##' @export
plot.riskIDM <- function(x, ...){
    require(ggplot2)
    out <- autoplot.riskIDM(x, ...)
    print(out)
    return(invisible(out))
}


## * autoplot.riskIDM (documentation)
##' @title Graphical display for For Illness Death Model
##' @description Diplay state occupancy probability of an Illness Death Model
##' @rdname autoplot.riskIDM
##'
##' @param object [riskIDM] output of the \code{riskIDM} function.
##' @param by [character] should the occupancy probabilities of all states be displayed on the same graphic, possibly using a different panel for each scenario (\code{"scenario"}).
##' Or should the occupancy probabilities for all scenarios be displayed on the same graphic, possibly using a different panel for each state (\code{"state"}).
##' @param scenario [character vector] name of the scenarios to be displayed. Use \code{"all"} to display all scenarios.
##' @param state [character vector] name of the states to be displayed. Use \code{"all"} to display all states.
##' @param indiv [logical or character vector] should the occupancy probabilities be displayed separately for each combination of covariates using a different type of line.
##' Not available for \code{stackplot=TRUE}.
##' @param ci [logical] should pointwise confidence intervals be displayed.
##' Not available for \code{stackplot=TRUE} and require a non-parametric bootstrap has been performed when running  the \code{riskIDM} function.
##' @param stackplot [logical] should the occupancy probability be cumulated over states under a specific scenario?
##' Only relevant when \code{by} equals \code{"scenario"}.
##' @param linewidth [numeric, >0] thickness of the line used to display the occupancy probabilities.
##' Only relevant for \code{type="curve"} and \code{type="stackcurve"}.
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
##' plot(e.riskPH)
##' plot(e.riskPH, scenario = "all")
##' plot(e.riskPH, stackplot = FALSE)
##' plot(e.riskPH2, stackplot = FALSE, scenario = "all")
##' 
##' plot(e.riskPH, by = "state")
##' plot(e.riskPH, by = "state", state = "all")
##' plot(e.riskPH2, by = "state", state = "all")
##'
##' plot(e.riskPH, by = "contrast")
##' 
##' #### fit IDM with covariates ####
##' e.riskPH3 <- riskIDM(~age, data = ebmt3, PH = FALSE, 
##'                      var.id = "id", keep.indiv = TRUE,
##'                      var.type = c("prstat", "rfsstat"),
##'                      var.time = c("prtime", "rfstime"))
##' 
##' plot(e.riskPH3, by = "state", state = "all", indiv = TRUE)
##' plot(e.riskPH3, by = "state", state = "all", indiv = ">40")
##' plot(e.riskPH3, by = "state", state = "all", indiv = c(">40","<=20"))
##' plot(e.riskPH3, by = "state", state = "all", indiv = FALSE)



## * autoplot.riskIDM (code)
##' @name autoplot.riskIDM
##' @export
autoplot.riskIDM <- function(object, by = "scenario", scenario = NULL, state = NULL, indiv = FALSE, ci = NULL,
                             stackplot = NULL, metric = "difference", linewidth = 2, ci.alpha = 0.2, breaks = NULL, ...){

    ## ** normalize arguments
    by <- match.arg(by, c("scenario","state","contrast"))
    
    dots <- list(...)
    if (length(dots) > 0) {
        stop("Unknown argument(s) '", paste(names(dots), collapse = "' '"), "'. \n")
    }

    n.boot <- object$args$n.boot
    if(!identical(indiv,FALSE)){
        if(!is.null(ci) && ci>0){
            stop("Incorrect argument \'ci\': cannot display uncertainty for individual occupancy probabilities. \n")
        }else{
            ci <- FALSE
        }
        if(object$args$indiv == FALSE){
            stop("Incorrect value for argument \'indiv\': individual occupancy probabilities have not been stored. \n",
                 "Consider setting the argument \'indiv\' to TRUE when calling riskIDM. \n")
        }
        if(!is.null(stackplot) && stackplot == TRUE){
            warning("Argument \'indiv\' ignored when argument \'stackplot\' is TRUE. \n")
            indiv <- FALSE
        }
        if(is.character(indiv)){
            indiv <- match.arg(indiv, unique(object$jump.indiv$signature), several.ok = TRUE)
        }
    }else if(!is.null(ci) && ci>0 && n.boot>0){
        stop("Incorrect argument \'ci\': require to have performed non-parametric bootstrap to display uncertainty. \n",
             "Consider setting the argument \'n.boot\' to a large value (e.g. 10000) when calling riskIDM. \n")
       
    }else{
        ci <- n.boot>0
    }
    if(is.null(stackplot)){
        stackplot <- (by == "scenario")
    }else if(stackplot && by == "state"){
        stop("Incorrect value for argument stackplot: should be FALSE when argument \'by\' is \"state\". \n")
    }
    object.scenario <- object$scenario
    if(is.null(scenario)){
        if(by == "scenario"){scenario <- "observed"}else{scenario <- object.scenario}
    }else if(identical(scenario,"all")){
        scenario <- object.scenario
    }else{
        scenario <- match.arg(scenario, object.scenario, several.ok = TRUE)
    }
    object.states <- object$states
    if(is.null(state)){
        if(by == "scenario"){state <- object.states$name}else{state <- utils::tail(object.states$name,1)}
    }else if(identical(state,"all")){
        state <- object.states$name
    }else{
        state <- match.arg(state, object.states$name)
    }
    if(is.null(breaks) && by!="contrast"){
        breaks  <- seq(0,1,0.1)
    }

    ## ** dataset
    if(!identical(indiv,FALSE)){
        table.gg <- data.table::copy(object$jump.indiv)
        if(is.character(indiv)){
            table.gg <- table.gg[table.gg$signature %in% indiv,]
        }
        indiv.levels <- unique(table.gg$signature)
        indiv.linetype <- setNames(1 + 1:length(indiv.levels), indiv.levels)
    }else if(by == "contrast"){
        table.gg <- model.tables(object, contrast = scenario, state = state, metric = metric)
        if(metric == "difference"){
            table.gg$scenario <- paste0(table.gg$alternative," - ",table.gg$reference)
        }else if(metric == "ratio"){
            table.gg$scenario <- paste0(table.gg$alternative," / ",table.gg$reference)
        }
        indiv.levels <- NULL
        indiv.linetype <- NULL
        scenario <- unique(table.gg$scenario)
    }else{
        table.gg <- data.table::copy(object$jump.estimate)
        indiv.levels <- NULL
        indiv.linetype <- NULL
    }
    ## subset in two steps to avoid confusion between the argument (state,scenario) and the column names (since data.table is used)
    index.keep <- intersect(which(table.gg$state %in% state),which(table.gg$scenario %in% scenario))
    table.gg <- table.gg[intersect(index.keep,which(!is.na(table.gg$estimate)))]
    table.gg$scenario <- factor(table.gg$scenario, scenario)

    if(stackplot){
        tol <- object$tol
        if(!identical(indiv,FALSE)){
            table.gg <- do.call(rbind,by(table.gg,interaction(table.gg$scenario,table.gg$state,table.gg$signature), function(iDF){
                iDF2 <- iDF[1:(NROW(iDF)-1)]
                iDF2$time <- iDF$time[2:NROW(iDF)]-tol
                iOut <- rbind(iDF,iDF2)
                setkeyv(iOut, c("scenario","state","signature","time"))
                return(iOut)
            }))
        }else{
            table.gg <- do.call(rbind,by(table.gg,interaction(table.gg$scenario,table.gg$state), function(iDF){
                iDF2 <- iDF[1:(NROW(iDF)-1)]
                iDF2$time <- iDF$time[2:NROW(iDF)]-tol
                iOut <- rbind(iDF,iDF2)
                setkeyv(iOut, c("scenario","state","time"))
                return(iOut)
            }))
        }
    }
    ## ** prepare for graph 
    Ustate <- unique(table.gg$state)
    Uscenario <- unique(table.gg$scenario)
    if(by == "scenario"){
        name.color <- c("state","State")
        name.facet <- "scenario"
        if(length(Uscenario)==1){
            label.y <- paste0("Occupancy probability for scenario \"",Uscenario,"\"")
            formula.facet <- NULL
        }else{
            label.y <- "Occupancy probability"
            formula.facet <- ~scenario
        }
        if(is.numeric(object.states$all)){
            table.gg$state <- factor(table.gg$state, levels = object.states$name, labels = object.states$all)
        }
    }else if(by == "state"){
        name.color <- c("scenario","Scenario")
        name.facet <- "state"
        if(length(Ustate)==1){
            if(is.numeric(object.states$all)){
                label.y <- paste0("Occupancy probability for state ",object.states$all[object.states$name==Ustate])
            }else{
                label.y <- paste0("Occupancy probability for state \"",Ustate,"\"")
            }
            formula.facet <- NULL
        }else{
            label.y <- "Occupancy probability"
            formula.facet <- ~state
        }
        if(is.numeric(object.states$all)){
            table.gg$state <- factor(table.gg$state, levels = object.states$name, labels = paste0("state ",object.states$all))
        }
    }else if(by == "contrast"){
        name.color <- c("scenario","Contrast")
        name.facet <- "state"
        if(is.numeric(object.states$all)){
            label.y <- paste0("Occupancy probability for state ",object.states$all[object.states$name==Ustate])
        }else{
            label.y <- paste0("Occupancy probability for state \"",Ustate,"\"")
        }
        formula.facet <- NULL
        if(is.numeric(object.states$all)){
            table.gg$state <- factor(table.gg$state, levels = object.states$name, labels = paste0("state ",object.states$all))
        }
    }
    if(stackplot){
        table.gg$state <- factor(table.gg$state, levels = rev(levels(table.gg$state)))
    }
    
    ## ** graphical display
    if(!identical(indiv,FALSE)){
        gg <- ggplot2::ggplot(table.gg, ggplot2::aes(x = time, y = estimate, linetype = signature, group = interaction(.data[[name.color[1]]],signature,drop=TRUE)))
        gg <- gg + ggplot2::scale_linetype_manual(values = indiv.linetype, breaks = names(indiv.linetype))            
    }else{
        gg <- ggplot2::ggplot(table.gg, ggplot2::aes(x = time, y = estimate, group = .data[[name.color[1]]]))
        if(stackplot == FALSE & n.boot>0 & ci){
            gg <- gg + ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, fill = .data[[name.color[1]]]), alpha = ci.alpha)
            gg <- gg + ggplot2::labs(fill = name.color[2])
        }
    }
    if(stackplot){
        gg <- gg + ggplot2::geom_area(ggplot2::aes(fill = .data[[name.color[1]]])) 
        gg <- gg + ggplot2::labs(fill = name.color[2], y = label.y)
        gg <- gg + ggplot2::coord_cartesian(ylim = c(0,1))
    }else{
        gg <- gg + ggplot2::geom_step(linewidth = linewidth, ggplot2::aes(y = estimate, color = .data[[name.color[1]]]))
        gg <- gg + ggplot2::labs(color = name.color[2], y = label.y)
    }
    if(is.null(breaks)){
        gg <- gg + ggplot2::scale_y_continuous(labels=scales::percent)
    }else{
        gg <- gg + ggplot2::scale_y_continuous(breaks = breaks, labels=scales::percent)
    }
    if(!is.null(formula.facet)){
        gg <- gg + ggplot2::facet_wrap(formula.facet)
    }
    gg <- gg + ggplot2::theme(axis.ticks.length=unit(.25, "cm"),
                              legend.key.width = unit(3,"line"))


    ## ** export
    return(gg)
}

## * .hazard2risk
.hazard2risk <- function(time, hazard01, hazard02, hazard12, states, prodlim = TRUE){

    n.time <- length(time)
    n.switch <- NCOL(hazard01)

    cumhazard01 <- colCumSum(hazard01)
    cumhazard02 <- cumsum(hazard02)
    cumhazard12 <- colCumSum(hazard12)

    ## ** pre-compute
    if(prodlim){
        S11 <- cumprod(1-rowSums(hazard01)-hazard02)
    }else{
        S11 <- exp(- rowSums(cumhazard01) - cumhazard02)
    }
    S11m <- c(1,S11[-n.time])
    S11m.hazard01 <- riskRegression::colMultiply_cpp(hazard01, S11m)

    if(prodlim){
        ## ## fast implementation
        C1.cumhazard12 <- apply(1-hazard12,2,cumprod)
        S11m.hazard01.C1.cumhazard12.inv <- riskRegression::colCumSum(S11m.hazard01/C1.cumhazard12)
        S01 <- S11m.hazard01.C1.cumhazard12.inv * C1.cumhazard12

        ## ## ## slow but explicit implementation
        ## range(S01 - do.call(rbind,lapply(1:n.time, function(iTau){ ## iTau <- 10
        ##     if(iTau==1){return(0)} ## both hazard cannot be simulataneously non-0 at time 1
        ##     iScale <- matrix(apply(1-hazard12[1:iTau,,drop=FALSE],2,prod), nrow = iTau, ncol = n.switch, byrow = TRUE)
        ##     iS22 <- iScale/apply(1-hazard12[1:iTau,,drop=FALSE],2,cumprod)            
        ##     iFactor <- matrix(S11m[1:iTau], nrow = iTau, ncol = n.switch, byrow = FALSE)
        ##     return(colSums(iFactor * hazard01[1:iTau,,drop=FALSE] * iS22[1:iTau,,drop=FALSE]))
        ## })))

    }else{ 
        ## ## fast implementation
        e.cumhazard12 <- exp(cumhazard12)
        S11m.hazard01.e.cumhazard12 <-  riskRegression::colCumSum(S11m.hazard01 * e.cumhazard12)
        S01 <- S11m.hazard01.e.cumhazard12 / e.cumhazard12

        ## ## slow but explicit implementation
        ## range(S01 - do.call(rbind,lapply(1:n.time, function(iTau){ ## iTau <- 10
        ##     if(iTau==1){return(0)} ## both hazard cannot be simulataneously non-0 at time 1
        ##     iCenter <- matrix(cumhazard12[iTau,,drop=FALSE], nrow = iTau, ncol = n.switch, byrow = TRUE)
        ##     iS22 <- exp(-iCenter+cumhazard12[1:iTau,,drop=FALSE])
        ##     iFactor <- matrix(S11m[1:iTau], nrow = iTau, ncol = n.switch, byrow = FALSE)
        ##     return(colSums(iFactor * hazard01[1:iTau,,drop=FALSE] * iS22[1:iTau,,drop=FALSE]))
        ## })))
    }
    S01m <- rbind(rep(0,n.switch),S01[-n.time,,drop=FALSE])
    
    ## ** scenario
    out <- data.frame(time,
                      S11,
                      S01,
                      cumsum(hazard02 * S11m) + rowSums(colCumSum(hazard12 * S01m))
                      )
    
    ## ** export
    names(out) <- c("time",states$name)
    return(out)
}

## * .commonString
## find the common consecutive part between two strings
## copied from https://stackoverflow.com/questions/35381180/identify-a-common-pattern
## .commonString("aaachang2","aaabbb")
## .commonString("aaa55change2","aaachangebbb")
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
