# {{{ doc
#' @title Survival curve using ggplot2
#' @description Display a non-parametric or semi-parametric survival curve. 
#' @name ggSurv
#' 
#' @param object a coxph object or a survfit object or a data.table object.
#' @param newdata dataset containing the covariates conditional to which the survival is displayed.
#' @param timeVar the name of the column in the data.table object containing the time variable
#' @param survivalVar the name of the column in the data.table object containing the survival at each time
#' @param ciInfVar the name of the column in the data.table object containing the lower bound of the confidence interval 
#' @param ciSupVar the name of the column in the data.table object containing the upper bound of the confidence interval 
#' @param eventVar the name of the column in the data.table object containing the number of events at each time
#' @param censorVar the name of the column in the data.table object containing the number of censored patients at each time
#' @param strataVar the name of the column in the data.table object containing the strata variable
#' @param format the format used to export the data. Can be data.table or data.frame.
#' @param plot should the plot be displayed
#' @param legend.position,title,textSize annotations on the plot
#' @param ylim y range (survival)
#' @param line.size size of the survival curves
#' @param confint,alpha.CIfill,alpha.CIline how to display confidence intervals
#' @param censoring,alpha.censoring,colour.censoring,shape.censoring,size.censoring,name.censoring how to display censored events
#' @param events,alpha.events,colour.events,shape.events,size.events,name.events how to display non-censored events
#' @param ... arguments to be passed to lower level functions
#'
#' @return A list:
#' \itemize{
#'   \item plot: the ggplot object
#'   \item data: the data used to create the ggplot object
#' }
#' @examples
#' library(lava)
#' library(survival)
#'
#' set.seed(10)
#' n <- 500
#' newdata <- data.frame(X1=1)
#' time <- 0.25
#'
#' ## simulate non proportional hazard using lava
#' m <- lvm()
#' regression(m) <- y ~ 1
#' regression(m) <- s ~ exp(-2*X1)
#' distribution(m,~X1) <- binomial.lvm()
#' distribution(m,~cens) <- coxWeibull.lvm(scale=1)
#' distribution(m,~y) <- coxWeibull.lvm(scale=1,shape=~s)
#' eventTime(m) <- eventtime ~ min(y=1,cens=0)
#' d <- as.data.table(sim(m,n))
#' setkey(d, eventtime)
#'
#' ## Cox (PH)
#' m.cox <- coxph(Surv(eventtime, status) ~ X1, data = d, y = TRUE, x = TRUE)
#' res1 <- ggSurv(m.cox)
#' ggSurv(m.cox, newdata = d[,.SD[1], by = X1])
#' ggSurv(m.cox, newdata = d[,.SD[1], by = X1], confint = TRUE)
#'
#' ## Cox (stratified)
#' mStrata.cox <- coxph(Surv(eventtime, status) ~ strata(X1), data = d, y = TRUE, x = TRUE)
#' ggSurv(mStrata.cox)
#'
#' ## KM
#' m.KM <- survfit(Surv(eventtime, status) ~ X1, data = d)
#' res2 <- ggSurv(m.KM)
#'
#' ## combine plot
#' dt.gg <- rbind(cbind(res1$data[original == TRUE & X1 == 1], model = "cox"),
#'               cbind(res2$data[original == TRUE & X1 == 1], model = "KM")
#' )
#' dt.gg[,c("X1","original") := NULL]
#' ggSurv(dt.gg, timeVar = "time", survivalVar = "survival", strataVar = "model")
#'
#' @export
`ggSurv` <-
  function(object,...) UseMethod("ggSurv")
# }}}

# {{{ ggSurv.survfit
#' @rdname ggSurv
#' @export
ggSurv.survfit <- function(object, ...){
  
  dt.ggplot <- data.table::data.table(time = object$time,
                                      survival = object$surv,
                                      ci.inf = object$lower,
                                      ci.sup = object$upper,
                                      n.event = object$n.event,
                                      n.censor = object$n.censor)
  
  strata <- object$strata
  
  ## create strata 
  if(!is.null(strata)){
    Flevels.strata <- names(strata)
    base.strata <- sapply(Flevels.strata, function(s){unlist(lapply(strsplit(s, split = "=", fixed = TRUE), "[",1))})
    varStrata <- unique(base.strata)
    
    levels.model <- mapply(FUN = function(base, i){gsub(pattern = paste0(base,"="), replacement = "", x = i, fixed = TRUE)}, 
                           base = base.strata, i = Flevels.strata)
    n.levels.model <- length(levels.model)
    
    valueStrata <-  unlist( lapply(1:n.levels.model, function(s){rep(levels.model[s], times = strata[s])}) )
    dt.ggplot[, (varStrata) := valueStrata]
  }else{
    varStrata <- NULL 
  }
    setkey(dt.ggplot, time)

    res <- ggSurv(object = dt.ggplot,
                  timeVar = "time", survivalVar = "survival", ciInfVar = "ci.inf", ciSupVar = "ci.sup", 
                  eventVar = "n.event", censorVar = "n.censor",  strataVar = varStrata,
                  ...)
    return(invisible(res))
}
# }}}

# {{{ ggSurv.coxph
#' @rdname ggSurv
#' @export
ggSurv.coxph <- function(object, newdata = NULL, confint = FALSE, ...){

    ## for CRAN test
    . <- Utimes <- xxSTRATAxx <- n.event <- n.censor <- survival.lower <- survival.upper <- strata <- status <- survival <- NULL

    ## find all event times
    data <- as.data.table(riskRegression::CoxDesign(object))[,.SD,.SDcols=c("stop","status")]

    ## find all strata in the original object
    coxInfo <- riskRegression::CoxVariableName(object)
    name.Xstrata <- coxInfo$stratavars.original
    #terms.special <- prodlim::strip.terms(terms(CoxFormula(object)), special = "strata")
    name.Xlp <- coxInfo$lpvars#attr(terms.special, "term.labels")
    name.X <- c(name.Xstrata,name.Xlp)

    ## initialize newdata
    if(is.null(newdata)){
        originalData <- as.data.table(eval(object$call$data))
        setnames(originalData, old = coxInfo$time, new = "stop")
        
        if(length(name.X)>0){
            newdata <- originalData[,.SD,.SDcols=c("stop",name.X)]
        }else{
            newdata <- originalData[,.SD,.SDcols = "stop"]
        }        
        init.newdata <- TRUE
    }else{
        init.newdata <- FALSE
    }

    ## add strata
    if(length(name.X)==0){
        newdata[,xxSTRATAxx := "1"]
        strataVar <- NULL
    }else if(length(name.X)==1){
        newdata[,xxSTRATAxx := as.character(.SD[[1]]),.SDcols = name.X]
        strataVar <- name.X
    }else{
        newdata[,xxSTRATAxx := interaction(.SD),.SDcols = name.X]
        strataVar <- paste(name.X, collpase=".", sep = "")
    }
    
    ## add prediction times
    if(init.newdata){ # also restrict to one observation per strata
        newdata <- newdata[,cbind(Utimes = .(.(unique(stop))),.SD[1]),by = xxSTRATAxx]        
    }else{
        newdata[,Utimes := .(.(sort(unique(data$stop))))]        
    }

    ## compute predictions
    n.newdata <- NROW(newdata)    
    predC <- NULL    
    for(iObs in 1:n.newdata){ # iObs <- 1
        rrPred  <- as.data.table(riskRegression::predictCox(object, newdata = newdata[iObs], times = newdata[iObs,Utimes[[1]]],
                                                            se = confint, keep.strata = FALSE, type = "survival"))
 
        predC <- rbind(predC,
                       cbind(rrPred, strata = newdata[iObs,xxSTRATAxx]) )
    }
    setnames(predC, old = c("times","strata"), new = c("time",strataVar))

    ## normalize for export
    predC[, n.event := NA]
    predC[, n.censor := NA]
    if(confint==FALSE){
        predC[, survival.lower := NA]
        predC[, survival.upper := NA]
    }

    ## export
    res <- ggSurv(object = predC,
                  timeVar = "time", survivalVar = "survival",
                  confint = confint,
                  ciInfVar = "survival.lower",
                  ciSupVar = "survival.upper",
                  eventVar = "n.event", censorVar = "n.censor",  strataVar = strataVar, 
                  ...)
  
    return(invisible(res))
  
}
# }}}

# {{{ ggSurv.data.table
#' @rdname ggSurv
#' @export
ggSurv.data.table <- function(object, format = "data.table",
                              timeVar = "time", survivalVar = "survival", ciInfVar = NULL, ciSupVar = NULL, 
                              eventVar = NULL, censorVar = NULL,  strataVar = NULL, 
                              plot = TRUE, legend.position = "top",
                              title = NULL, textSize = NULL, ylim = NULL, line.size = 2,
                              confint = FALSE, alpha.CIfill = 0.2, alpha.CIline = 0.5,
                              censoring = FALSE, alpha.censoring = 1, colour.censoring = rgb(0.2,0.2,0.2), shape.censoring = 3, size.censoring = 2, name.censoring = "censoring",
                              events = FALSE, alpha.events = 1, colour.events = rgb(0.2,0.2,0.2), shape.events = 8, size.events = 2, name.events = "event", ...){

    ## for CRAN test
    original <- n.censor <- n.event <- survival <- ci.inf <- ci.sup <- NULL
    
    ## avoid unexpected modification of the original object
    object <- copy(object)

    ## remove useless columns in object
    object.name <- names(object)
    if(confint==FALSE && !is.null(ciInfVar) && ciInfVar %in% object.name){
        object[,(ciInfVar) := NULL]
        ciInfVar <- NULL
    }
    if(confint==FALSE && !is.null(ciSupVar) && ciSupVar %in% object.name){
        object[,(ciSupVar) := NULL]
        ciSupVar <- NULL
    }
    if(censoring==FALSE && !is.null(censorVar) && censorVar %in% object.name){
        object[,(censorVar) := NULL]
        censorVar <- NULL
    }
    if(events==FALSE && !is.null(eventVar) && eventVar %in% object.name){
        object[,(eventVar) := NULL]
        eventVar <- NULL
    }

    ## check object
    butils.base::validCharacter(value1 = timeVar, name1 = "timeVar", validLength = 1, validValues = object.name, method = "ggSurv.dt")
    butils.base::validCharacter(value1 = survivalVar, name1 = "survivalVar", validLength = 1, validValues = object.name, method = "ggSurv.dt")
    butils.base::validCharacter(value1 = ciInfVar, name1 = "ciInfVar", validLength = 1, refuse.NULL = confint, validValues = object.name, method = "ggSurv.dt")
    butils.base::validCharacter(value1 = ciSupVar, name1 = "ciSupVar", validLength = 1, refuse.NULL = confint, validValues = object.name, method = "ggSurv.dt")
    butils.base::validCharacter(value1 = eventVar, name1 = "eventVar", validLength = 1, refuse.NULL = events, validValues = object.name, method = "ggSurv.dt")
    butils.base::validCharacter(value1 = censorVar, name1 = "censorVar", validLength = 1, refuse.NULL = censoring, validValues = object.name, method = "ggSurv.dt")
    butils.base::validCharacter(value1 = strataVar, name1 = "strataVar", validLength = 1, refuse.NULL = FALSE, validValues = object.name, method = "ggSurv.dt")

    ## normalize object
    object <- object[,.SD, .SDcols = c(timeVar,survivalVar,ciInfVar,ciSupVar, eventVar, censorVar, strataVar)]
    setnames(object, old = c(timeVar,survivalVar), new = c("time","survival"))

    if(confint){
        setnames(object, old = c(ciInfVar,ciSupVar), new = c("ci.inf","ci.sup"))        
    }else{
        object[,ci.inf := NA]
        object[,ci.sup := NA]
    }
    ciInfVar <- "ci.inf"
    ciSupVar <- "ci.sup"
    
    if(censoring){
        setnames(object, old = eventVar, new = "n.event")
    }else{
        object[,n.event := NA]
    }
    eventVar <- "n.event"
    
    if(events){
        setnames(object, old = censorVar, new = "n.censor")
    }else{
        object[,n.censor := NA]
    }
    censorVar <- "n.censor"
    
    ## order dataset
    setkeyv(object, c("time", strataVar))

    #### duplicate dataset to add
    ## initial point  (t=0)
    ## point just before the jump time
    dt0 <- data.table(survival = 1, ci.inf = NA, ci.sup = NA)

    if (is.null(strataVar)) {
        copyx <- rbind(dt0,
                       object[-.N,.SD[, c(survivalVar, ciInfVar, ciSupVar), with = FALSE]]
                       )
        copyx[, time := object[,time - .Machine$double.eps*100]]
        copyx <- rbind(cbind(time = 0, dt0), copyx)
    }else{
        levels.model <- unique(object[[strataVar]])
        n.levels.model <- length(levels.model)
        copyx <- NULL
        
        for(iStrata in 1:n.levels.model){ # iStrata <- 1
            iLevel <- levels.model[iStrata]
            iIndex <- which(object[[strataVar]]==iLevel)
            iDt0 <- copy(dt0) ; iDt0[,(strataVar) := iLevel] ; 

            iCopyx <- rbind(iDt0,
                            object[iIndex[-length(iIndex)],.SD[, c(survivalVar, ciInfVar, ciSupVar, strataVar), with = FALSE]]
                            )
            iCopyx[, time := object[iIndex,time - .Machine$double.eps*100]]

            iCopyx <- rbind(cbind(iDt0, time = 0), iCopyx)
            copyx <- rbind(copyx, iCopyx)
        }
    }
    
    copyx[,n.censor := 0]
    copyx[,n.event := 0]
    object[,original := TRUE]
    copyx[,original := FALSE]

    setcolorder(copyx, names(object))   
    object <- rbind(object, copyx)
  
    #### basic display  
    if(!is.null(plot)){
        dt2.ggplot <- copy(object)
    
        if (is.null(strataVar)) {
            gg.base <- ggplot(data = dt2.ggplot, mapping = aes(x = time, y = survival))
        }else{
            gg.base <- ggplot(data = dt2.ggplot, mapping = aes_string(x = "time", y = "survival", color = strataVar, fill = strataVar))
        }
    
        # ci
        if(confint){
            if (is.null(strataVar)) {
                gg.base <- gg.base + geom_ribbon(aes(x = time, ymin = ci.inf, ymax = ci.sup, alpha = alpha.CIfill))
            }else{
                gg.base <- gg.base + geom_ribbon(data = dt2.ggplot, aes(x = time, ymin = ci.inf, ymax = ci.sup), alpha = alpha.CIfill)
                gg.base <- gg.base + geom_line(data = dt2.ggplot, aes(x = time, y = ci.inf), alpha = alpha.CIline)
                gg.base <- gg.base + geom_line(data = dt2.ggplot, aes(x = time, y = ci.sup), alpha = alpha.CIline)
            }
        }
    
        # step function
        gg.base <- gg.base + geom_line(data = dt2.ggplot, size = line.size)
        if (is.null(strataVar)) {
            gg.base <- gg.base + guides(fill = guide_legend(title = strataVar, title.position = legend.position),
                                        color = guide_legend(title = strataVar, title.position = legend.position)) #+ scale_colour_discrete(name = varStrata)
        }
    
        # censoring
        values.ShapeLegend <- vector(mode = "numeric", length =  events + censoring)
        names(values.ShapeLegend)  <- c(if(censoring){name.censoring}, if(events){"event"})
    
        if (censoring == TRUE) {
            eval(parse(text = paste0(
                           "gg.base <- gg.base + geom_point(data = dt2.ggplot[n.censor > 0], aes(x = time, y = survival, shape = \"",name.censoring,"\"), 
        alpha = alpha.censoring, ",if(!is.null(colour.censoring)){"colour = colour.censoring,"}," size = size.censoring)"
        )))
            values.ShapeLegend[name.censoring] <- shape.censoring
        } 
        # events
        if (events == TRUE) {
            eval(parse(text = paste0(
                           "gg.base <- gg.base + geom_point(data = dt2.ggplot[n.event > 0], aes(x = time, y = survival, shape = \"",name.events,"\"), 
        alpha = alpha.events, ",if(!is.null(colour.events)){"colour = colour.events,"}," size = size.events)"
        )))
            values.ShapeLegend[name.events] <- shape.events
        }
    
        if (censoring == TRUE || events == TRUE) {
            gg.base <- gg.base + scale_shape_manual( values = values.ShapeLegend )
            gg.base <- gg.base + guides(shape = guide_legend(title = "Observation", title.position = legend.position ))
        }
    
        if(!is.null(title)){gg.base <- gg.base + ggtitle(title)}
        if(!is.null(textSize)){gg.base <- gg.base +  theme(text = element_text(size = textSize))}
        if(!is.null(ylim)){gg.base <- gg.base + coord_cartesian(ylim = ylim)}
    
        if(plot==TRUE){
            print(gg.base)
        }
    }else{
        gg.base <- NULL
    }
  
    #### export
    if(format == "data.frame"){
        object <- as.data.frame(object)
    }
  
    return(invisible(list(data = object,
                          plot = gg.base)))
}
# }}}




