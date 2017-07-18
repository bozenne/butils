# {{{ doc
#' @title Survival curve using ggplot2
#' @description Display a non-parametric or semi-parametric survival curve. Can handle one strata variable.
#' @name ggSurv
#' 
#' @param x a coxph object or a survfit object or a data.table object.
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
#' 
#' @examples
#' library(survival)
#' dt <- as.data.table(aml)
#' dt[,x:=as.factor(x)]
#' 
#' ## survfit
#' KM <- survfit(Surv(time, status) ~ x, data = dt)
#' ggSurv(KM)
#' ggSurv(KM, confint = TRUE)
#' 
#' ## coxph
#' Cox <- coxph(Surv(time, status) ~ x, data = dt, x = TRUE)
#' ggSurv(Cox)
#' ggSurv(Cox, confint = TRUE)
#' 
#' ## data.table
#' dt2 <- data.table(time = 1:10, 
#'                   survival = seq(1, by = -0.01, length.out = 10), 
#'                   n.censor = 0, n.event = 1)
#' ggSurv(dt2)
#'
#' dt3 <- data.table(time = 1:10, 
#'                   survival = seq(1, by = -0.01, length.out = 10)
#' )
#' ggSurv.data.table(dt3)
#' 
#' @export
`ggSurv` <-
  function(x,...) UseMethod("ggSurv")
# }}}

# {{{ ggSurv.survfit
#' @rdname ggSurv
#' @export
ggSurv.survfit <- function(x, ...){
  
  dt.ggplot <- data.table::data.table(time = x$time,
                                      survival = x$surv,
                                      ci.inf = x$lower,
                                      ci.sup = x$upper,
                                      # std.err = x$std.err,
                                      # n.risk = x$n.risk,
                                      n.event = x$n.event,
                                      n.censor = x$n.censor)
  
  strata <- x$strata
  
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
  
  res <- ggSurv(x = dt.ggplot,
                timeVar = "time", survivalVar = "survival", ciInfVar = "ci.inf", ciSupVar = "ci.sup", 
                eventVar = "n.event", censorVar = "n.censor",  strataVar = varStrata,
                ...)
  return(invisible(res))
}
# }}}

# {{{ ggSurv.coxph
#' @rdname ggSurv
#' @export
ggSurv.coxph <- function(x, confint = FALSE, ...){

    ## for CRAN test
    strata <- status <- survival <- NULL

    ## get Design matrix
    data <- as.data.table(riskRegression::CoxDesign(x))[,.SD,.SDcols=c("stop","status")]

    ## add covariates in the original form
    coxInfo <- riskRegression::CoxVariableName(x)
    originalData <- as.data.table(eval(x$call$data))
        
    name.Xstrata <- coxInfo$stratavars.original
    if(length(name.Xstrata)>0){
        data <- cbind(data, originalData[,.SD,.SDcols=name.Xstrata])
    }
    terms.special <- prodlim::strip.terms(terms(CoxFormula(x)), special = "strata")
    name.Xlp <- attr(terms.special, "term.labels")
    if(length(name.Xlp)>0){
        data  <- cbind(data,originalData[,.SD,.SDcols=name.Xlp])
    }
    
    # order rows before extraction of non-duplicated observations
    # such that event will be extracted in priority 
    data[, status := 1-status]
    setkeyv(data, c("stop","status"))
    data[, status := 1-status]

    # remove duplicates
    data0 <- copy(data)
    data0[, status := NULL]
    index.duplicated <- which(duplicated(data0))
    index.Nduplicated <- setdiff(1:NROW(data), index.duplicated)

 
    data <- data[index.Nduplicated]
    data[, observation := 1:.N]
    time <- unique(data$stop)
    predC <- NULL
    
    if((length(name.Xstrata)+length(name.Xlp))>0){
        strataVar <- c(name.Xstrata, name.Xlp)
        X.levels <- data[,unique(.SD), .SDcols = strataVar]
        n.levels <- NROW(X.levels)

        xXXx <- x # avoid confusion in name
        rrPred <- data[, .(dtres = .(
                               cbind(strata = interaction(.SD[1,strataVar,with = FALSE]),
                                     as.data.table(riskRegression::predictCox(xXXx, newdata = .SD, times = time, se = confint, type = "survival")))
                           )),
                       by = strataVar, .SDcols = names(data)]
        predC <- rbindlist(rrPred[,dtres])

        if(length(strataVar)>1){
            strataVar <- paste(strataVar, collpase=".", sep = "")
        }
        setnames(predC, old = "strata", new = strataVar)
    }else{
        predC <- riskRegression::predictCox(x, newdata = data, times = time, se = confint, type = "survival")
        strataVar <- NULL
    }
    predC <- merge(predC, data[,.SD,.SDcols = c("status",'observation') ], by = "observation")
    predC[,observation := NULL]
    setnames(predC, old = "times", new = "time")
    
    #### reduce to unique time
    if(confint==FALSE){
        predC[, survival.lower := NA]
        predC[, survival.upper := NA]
    }
    predC <- predC[,list(survival = survival[1],
                         survival.upper = survival.upper[1],
                         survival.lower = survival.lower[1],
                         n.event = sum(status==1), n.censor = sum(status==0)), by = c("time",strataVar)]

    res <- ggSurv(x = predC,
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
ggSurv.data.table <- function(x, format = "data.table",
                              timeVar = "time", survivalVar = "survival", ciInfVar = NULL, ciSupVar = NULL, 
                              eventVar = NULL, censorVar = NULL,  strataVar = NULL, 
                              plot = TRUE, legend.position = "top",
                              title = NULL, textSize = NULL, ylim = NULL, line.size = 2,
                              confint = FALSE, alpha.CIfill = 0.2, alpha.CIline = 0.5,
                              censoring = FALSE, alpha.censoring = 1, colour.censoring = rgb(0.2,0.2,0.2), shape.censoring = 3, size.censoring = 2, name.censoring = "censoring",
                              events = FALSE, alpha.events = 1, colour.events = rgb(0.2,0.2,0.2), shape.events = 8, size.events = 2, name.events = "event", ...){

  ## for CRAN test
  original <- n.censor <- n.event <- survival <- ci.inf <- ci.sup <- NULL
  ##
  x <- copy(x)

  #### names
  butils.base::validCharacter(value1 = timeVar, name1 = "timeVar", validLength = 1, validValues = names(x), method = "ggSurv.dt")
  butils.base::validCharacter(value1 = survivalVar, name1 = "survivalVar", validLength = 1, validValues = names(x), method = "ggSurv.dt")
  butils.base::validCharacter(value1 = ciInfVar, name1 = "ciInfVar", validLength = 1, refuse.NULL = confint, validValues = names(x), method = "ggSurv.dt")
  butils.base::validCharacter(value1 = ciSupVar, name1 = "ciSupVar", validLength = 1, refuse.NULL = confint, validValues = names(x), method = "ggSurv.dt")
  butils.base::validCharacter(value1 = eventVar, name1 = "eventVar", validLength = 1, refuse.NULL = events, validValues = names(x), method = "ggSurv.dt")
  butils.base::validCharacter(value1 = censorVar, name1 = "censorVar", validLength = 1, refuse.NULL = censoring, validValues = names(x), method = "ggSurv.dt")
  butils.base::validCharacter(value1 = strataVar, name1 = "strataVar", validLength = 1, refuse.NULL = FALSE, validValues = names(x), method = "ggSurv.dt")

  x <- x[,.SD, .SDcols = c(timeVar,survivalVar,ciInfVar,ciSupVar, eventVar, censorVar, strataVar)]
  setnames(x, old = c(timeVar,survivalVar), new = c("time","survival"))

  if(confint){
      setnames(x, old = c(ciInfVar,ciSupVar), new = c("ci.inf","ci.sup"))
      ciInfVar <- "ci.inf"
      ciSupVar <- "ci.sup"
  }
  if(!is.null(eventVar)){
    setnames(x, old = eventVar, new = "n.event")
  }
  if(!is.null(censorVar)){
    setnames(x, old = censorVar, new = "n.censor")
  }

  #### order dataset
  setkeyv(x, c("time", strataVar))

    #### add first and last point
    x[, original := TRUE]
    if (is.null(strataVar)) {
        copyx <- x[, cbind(time = c(0,.SD$time[-.N] + .Machine$double.eps*100),
                           .SD[, c(survivalVar, ciInfVar, ciSupVar, eventVar, censorVar), with = FALSE],
                           original = FALSE)]
    
    }else{
        n.levels.model <- length(unique(x[[strataVar]]))
        x[, (strataVar) := lapply(.SD, as.factor), .SDcols = strataVar]
    
    copyx <- x[, cbind(time = c(0,.SD$time[-.N] + .Machine$double.eps*100),
                       .SD[, c(survivalVar, ciInfVar, ciSupVar, eventVar, censorVar), with = FALSE],
                       original = FALSE),
               by = strataVar]
  }
    if(censoring){
        copyx[,n.censor := 0]
    }
    if(events){
        copyx[,n.event := 0]
    }
    setcolorder(copyx, names(x))   
    x <- rbind(x, copyx)
  
  #### basic display  
  if(!is.null(plot)){
    dt2.ggplot <- copy(x)
    
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
    x <- as.data.frame(x)
  }
  
  return(invisible(list(data.ggplot = x,
                        ggplot = gg.base)))
}
# }}}




