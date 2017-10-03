# {{{ doc
#' @title Survival curve using ggplot2
#' @description Display a non-parametric or semi-parametric survival curve. 
#' @name ggSurv
#' 
#' @param object a coxph object or a survfit object or a data.table object.
#' @param data the data that have been used to fit the model.
#' @param newdata dataset containing the covariates conditional to which the survival is displayed.
#' @param var.time the name of the column in the data.table object containing the time variable
#' @param var.survival the name of the column in the data.table object containing the survival at each time
#' @param var.ci.inf the name of the column in the data.table object containing the lower bound of the confidence interval 
#' @param var.ci.sup the name of the column in the data.table object containing the upper bound of the confidence interval 
#' @param var.event the name of the column in the data.table object containing the number of events at each time
#' @param var.censor the name of the column in the data.table object containing the number of censored patients at each time
#' @param var.strata the name of the column in the data.table object containing the strata variable
#' @param format the format used to export the data. Can be data.table or data.frame.
#' @param plot should the plot be displayed
#' @param legend.position,title,text.size annotations on the plot
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
#' n <- 100
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
#' d <- as.data.table(lava::sim(m,n))
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
#' ggSurv(mStrata.cox, censoring = TRUE, event = TRUE)
#'
#' ## KM
#' m.KM <- survfit(Surv(eventtime, status) ~ X1, data = d)
#' res2 <- ggSurv(m.KM)
#'
#' ## combine plot
#' dt.gg <- rbind(cbind(res1$data[original == TRUE & X1 == 1], model = "cox"),
#'                cbind(res2$data[original == TRUE & X1 == 1], model = "KM")
#' )
#' dt.gg[,c("X1","original") := NULL]
#' ggSurv(dt.gg, var.time = "time", var.survival = "survival", var.strata = "model")
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
                var.time = "time", var.survival = "survival", var.ci.inf = "ci.inf", var.ci.sup = "ci.sup", 
                var.event = "n.event", var.censor = "n.censor",  var.strata = varStrata,
                ...)
  return(invisible(res))
}
# }}}

# {{{ ggSurv.coxph
#' @rdname ggSurv
#' @export
ggSurv.coxph <- function(object, data = NULL, newdata = NULL, confint = FALSE, ...){
  requireNamespace("riskRegression")
  
  ## for CRAN test
  . <- Utimes <- xxSTRATAxx <- n.event <- n.censor <- survival.lower <- survival.upper <- strata <- status <- survival <- NULL

  ## find all strata in the original object
  coxInfo <- riskRegression::coxVariableName(object)
  name.Xstrata <- coxInfo$strata.vars.original
  #terms.special <- prodlim::strip.terms(terms(coxFormula(object)), special = "strata")
  name.Xlp <- coxInfo$lpvars#attr(terms.special, "term.labels")
  name.X <- c(name.Xstrata,name.Xlp)
  
  ## extract the data used to fit the model
  if(is.null(data)){
    originalData <- extractData(object, force = TRUE, convert2dt = TRUE)
  }else{
    originalData <- copy(as.data.table(data))
    setnames(originalData, old = coxInfo$time, new = "stop")
  }
  
  ## initialize newdata
  if(is.null(newdata)){
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
    originalData[,xxSTRATAxx := "1"]
    var.strata <- NULL
  }else if(length(name.X)==1){
    newdata[,xxSTRATAxx := as.character(.SD[[1]]),.SDcols = name.X]
    originalData[,xxSTRATAxx := as.character(.SD[[1]]),.SDcols = name.X]
    var.strata <- name.X
  }else{
    newdata[,xxSTRATAxx := interaction(.SD),.SDcols = name.X]
    originalData[,xxSTRATAxx := interaction(.SD),.SDcols = name.X]
    var.strata <- paste(name.X, collpase=".", sep = "")
  }
  
  ## add prediction times
  if(init.newdata){ # also restrict to one observation per strata
    newdata <- newdata[,cbind(Utimes = .(.(unique(stop))),.SD[1]),by = xxSTRATAxx]        
  }else{
    newdata[,Utimes := .(.(sort(unique(originalData$stop))))]        
  }

  ## compute predictions
  n.newdata <- NROW(newdata)    
  predC <- NULL    
  for(iObs in 1:n.newdata){ # iObs <- 1
    
    rrPred  <- as.data.table(riskRegression::predictCox(object, newdata = newdata[iObs], times = newdata[iObs,Utimes[[1]]],
                                                        se = confint, keep.strata = FALSE, type = "survival"))
    
    dt.tempo <- cbind(rrPred, strata = newdata[iObs,xxSTRATAxx])
    
    # add the number of events
    dt.tempo2 <- originalData[xxSTRATAxx == newdata[iObs,xxSTRATAxx], .(n.censor = sum(status == 0), n.event = sum(status == 1)), by = stop]
    index.missing <- which(dt.tempo$times %in% dt.tempo2$stop == FALSE)
    if(length(index.missing)>0){
      dt.tempo2 <- rbind(dt.tempo2, cbind(stop = dt.tempo$times[index.missing], n.censor = 0, n.event = 0))
    }
    dt.tempo <- merge(x = dt.tempo, y = dt.tempo2, by.x = "times", by.y = "stop")
    predC <- rbind(predC, dt.tempo)
  }
  setnames(predC, old = c("times","strata"), new = c("time",var.strata))

  ## normalize for export
  if(confint==FALSE){
    predC[, survival.lower := NA]
    predC[, survival.upper := NA]
  }
  
  ## export
  res <- ggSurv(object = predC,
                var.time = "time", var.survival = "survival",
                confint = confint,
                var.ci.inf = "survival.lower",
                var.ci.sup = "survival.upper",
                var.event = "n.event", var.censor = "n.censor",  var.strata = var.strata, 
                ...)
  
  return(invisible(res))
  
}
# }}}

# {{{ ggSurv.data.table
#' @rdname ggSurv
#' @export
ggSurv.data.table <- function(object, format = "data.table",
                              var.time = "time", var.survival = "survival", var.ci.inf = NULL, var.ci.sup = NULL, 
                              var.event = NULL, var.censor = NULL,  var.strata = NULL, 
                              plot = TRUE, legend.position = "top",
                              title = NULL, text.size = NULL, ylim = NULL, line.size = 2,
                              confint = FALSE, alpha.CIfill = 0.2, alpha.CIline = 0.5,
                              censoring = FALSE, alpha.censoring = 1, colour.censoring = rgb(0.2,0.2,0.2), shape.censoring = 3, size.censoring = 2, name.censoring = "censoring",
                              events = FALSE, alpha.events = 1, colour.events = rgb(0.2,0.2,0.2), shape.events = 8, size.events = 2, name.events = "event", ...){
  
  ## for CRAN test
  original <- n.censor <- n.event <- survival <- ci.inf <- ci.sup <- NULL
  
  ## avoid unexpected modification of the original object
  object <- copy(object)
  
  ## remove useless columns in object
  object.name <- names(object)
  if(confint==FALSE && !is.null(var.ci.inf) && var.ci.inf %in% object.name){
    object[,(var.ci.inf) := NULL]
    var.ci.inf <- NULL
  }
  if(confint==FALSE && !is.null(var.ci.sup) && var.ci.sup %in% object.name){
    object[,(var.ci.sup) := NULL]
    var.ci.sup <- NULL
  }
  if(censoring==FALSE && !is.null(var.censor) && var.censor %in% object.name){
    object[,(var.censor) := NULL]
    var.censor <- NULL
  }
  if(events==FALSE && !is.null(var.event) && var.event %in% object.name){
    object[,(var.event) := NULL]
    var.event <- NULL
  }
  
  ## check object
  validCharacter(value1 = var.time, name1 = "var.time", valid.length = 1, valid.values = object.name, method = "ggSurv.dt")
  validCharacter(value1 = var.survival, name1 = "var.survival", valid.length = 1, valid.values = object.name, method = "ggSurv.dt")
  validCharacter(value1 = var.ci.inf, name1 = "var.ci.inf", valid.length = 1, refuse.NULL = confint, valid.values = object.name, method = "ggSurv.dt")
  validCharacter(value1 = var.ci.sup, name1 = "var.ci.sup", valid.length = 1, refuse.NULL = confint, valid.values = object.name, method = "ggSurv.dt")
  validCharacter(value1 = var.event, name1 = "var.event", valid.length = 1, refuse.NULL = events, valid.values = object.name, method = "ggSurv.dt")
  validCharacter(value1 = var.censor, name1 = "var.censor", valid.length = 1, refuse.NULL = censoring, valid.values = object.name, method = "ggSurv.dt")
  validCharacter(value1 = var.strata, name1 = "var.strata", valid.length = 1, refuse.NULL = FALSE, valid.values = object.name, method = "ggSurv.dt")
  
  ## normalize object
  object <- object[,.SD, .SDcols = c(var.time,var.survival,var.ci.inf,var.ci.sup, var.event, var.censor, var.strata)]
  setnames(object, old = c(var.time,var.survival), new = c("time","survival"))
  
  if(confint){
    setnames(object, old = c(var.ci.inf,var.ci.sup), new = c("ci.inf","ci.sup"))        
  }else{
    object[,ci.inf := NA]
    object[,ci.sup := NA]
  }
  var.ci.inf <- "ci.inf"
  var.ci.sup <- "ci.sup"
  
  if(censoring){
    setnames(object, old = var.censor, new = "n.censor")
  }else{
    object[,n.censor := NA]
  }
  var.censor <- "n.censor"
  
  if(events){
    setnames(object, old = var.event, new = "n.event")
  }else{
    object[,n.event := NA]
  }
  var.event <- "n.event"
  
  ## order dataset
  setkeyv(object, c("time", var.strata))
  
  #### duplicate dataset to add
  ## initial point  (t=0)
  ## point just before the jump time
  dt0 <- data.table(survival = 1, ci.inf = NA, ci.sup = NA)
  
  if (is.null(var.strata)) {
    copyx <- rbind(dt0,
                   object[-.N,.SD[, c(var.survival, var.ci.inf, var.ci.sup), with = FALSE]]
    )
    copyx[, time := object[,time - .Machine$double.eps*100]]
    copyx <- rbind(cbind(time = 0, dt0), copyx)
  }else{
    levels.model <- unique(object[[var.strata]])
    n.levels.model <- length(levels.model)
    copyx <- NULL
    
    for(iStrata in 1:n.levels.model){ # iStrata <- 1
      iLevel <- levels.model[iStrata]
      iIndex <- which(object[[var.strata]]==iLevel)
      iDt0 <- copy(dt0) ; iDt0[,(var.strata) := iLevel] ; 
      
      iCopyx <- rbind(iDt0,
                      object[iIndex[-length(iIndex)],.SD[, c(var.survival, var.ci.inf, var.ci.sup, var.strata), with = FALSE]]
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
    
    if (is.null(var.strata)) {
      gg.base <- ggplot(data = dt2.ggplot, mapping = aes(x = time, y = survival))
    }else{
      gg.base <- ggplot(data = dt2.ggplot, mapping = aes_string(x = "time", y = "survival", color = var.strata, fill = var.strata))
    }
    
    # ci
    if(confint){
      if (is.null(var.strata)) {
        gg.base <- gg.base + geom_ribbon(aes(x = time, ymin = ci.inf, ymax = ci.sup, alpha = alpha.CIfill))
      }else{
        gg.base <- gg.base + geom_ribbon(data = dt2.ggplot, aes(x = time, ymin = ci.inf, ymax = ci.sup), alpha = alpha.CIfill)
        gg.base <- gg.base + geom_line(data = dt2.ggplot, aes(x = time, y = ci.inf), alpha = alpha.CIline)
        gg.base <- gg.base + geom_line(data = dt2.ggplot, aes(x = time, y = ci.sup), alpha = alpha.CIline)
      }
    }
    
    # step function
    gg.base <- gg.base + geom_line(data = dt2.ggplot, size = line.size)
    if (is.null(var.strata)) {
      gg.base <- gg.base + guides(fill = guide_legend(title = var.strata, title.position = legend.position),
                                  color = guide_legend(title = var.strata, title.position = legend.position)) #+ scale_colour_discrete(name = varStrata)
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
    if(!is.null(text.size)){gg.base <- gg.base +  theme(text = element_text(size = text.size))}
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




