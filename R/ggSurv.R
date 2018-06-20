## * ggSurv (documentation)
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
#' @param name.strata the title for the color and fill legend in the plot.
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
#'

## * ggSurv (examples)
#' @rdname ggSurv
#' @examples
#' library(lava)
#' library(survival)
#' library(data.table)
#' 
#' set.seed(10)
#' n <- 100
#' newdata <- data.frame(X1=1)
#' time <- 0.25
#'
#' ## simulate non proportional hazard using lava
#' m <- lvm()
#' regression(m) <- y ~ 1
#' regression(m) <- s ~ exp(-2*X1+X2)
#' distribution(m,~X1+X2) <- binomial.lvm()
#' distribution(m,~cens) <- coxWeibull.lvm(scale=1)
#' distribution(m,~y) <- coxWeibull.lvm(scale=1,shape=~s)
#' eventTime(m) <- eventtime ~ min(y=1,cens=0)
#' d <- as.data.table(lava::sim(m,n))
#' setkey(d, eventtime)
#'
#' ## Cox (PH)
#' ## no regressor
#' m.cox <- coxph(Surv(eventtime, status) ~ 1, data = d, y = TRUE, x = TRUE)
#' res1 <- ggSurv(m.cox)
#'
#' ## regressors
#' m.cox <- coxph(Surv(eventtime, status) ~ X1, data = d, y = TRUE, x = TRUE)
#' res1 <- ggSurv(m.cox)
#' ggSurv(m.cox, newdata = d[,.SD[1], by = X1])
#' ggSurv(m.cox, newdata = d[,.SD[1], by = X1], confint = TRUE)
#'
#' ## Cox (stratified)
#' mStrata.cox <- coxph(Surv(eventtime, status) ~ strata(X1), data = d, y = TRUE, x = TRUE)
#' ggSurv(mStrata.cox, censoring = TRUE, event = TRUE)
#' 
#' mStrata.cox <- coxph(Surv(eventtime, status) ~ strata(X1) + strata(X2), data = d, y = TRUE, x = TRUE)
#' ggSurv(mStrata.cox, censoring = TRUE, event = TRUE)
#' 
#' ## KM
#' m.KM <- survfit(Surv(eventtime, status) ~ X1, data = d)
#' res2 <- ggSurv(m.KM)
#'
#' ## combine plot
#' res1$data[,X1 := as.numeric(factor(strata2, levels = c("X1=0","X1=1")))-1]
#' res1$data[,strata2 := NULL]
#' dt.gg <- rbind(cbind(res1$data[original == TRUE & X1==1], model = "cox"),
#'                cbind(res2$data[original == TRUE & X1==1], model = "KM")
#' )
#' dt.gg[,c("X1","original") := NULL]
#' ggSurv(dt.gg, var.time = "time", var.survival = "survival", var.strata = "model")
#'
#' @export
`ggSurv` <-
  function(object,...) UseMethod("ggSurv")

## * ggSurv.survfit (code)
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
    nameStrata <- "strata"
    
    levels.model <- mapply(FUN = function(base, i){gsub(pattern = paste0(base,"="), replacement = "", x = i, fixed = TRUE)}, 
                           base = base.strata, i = Flevels.strata)
    n.levels.model <- length(levels.model)
    
    valueStrata <-  unlist( lapply(1:n.levels.model, function(s){rep(levels.model[s], times = strata[s])}) )
    dt.ggplot[, (varStrata) := valueStrata]
  }else{
    varStrata <- NULL 
    nameStrata <- NULL 
  }
  setkey(dt.ggplot, time)
  
  res <- ggSurv(object = dt.ggplot,
                var.time = "time", var.survival = "survival", var.ci.inf = "ci.inf", var.ci.sup = "ci.sup", 
                var.event = "n.event", var.censor = "n.censor",  name.strata = nameStrata, var.strata = varStrata,
                ...)
  return(invisible(res))
}


## * ggSurv.coxph (code)
#' @rdname ggSurv
#' @export
ggSurv.coxph <- function(object, data = NULL, newdata = NULL, confint = FALSE, ...){

    tryPkg <- requireNamespace("riskRegression")
    if("try-error" %in% class(tryPkg)){
        stop(tryPkg)
    }
  
    ## for CRAN test
    . <- strata <- NULL

    ## extract the data used to fit the model
    object.data <- extractData(object, design.matrix = TRUE)

    ## find all strata in the original object
    coxInfo <- riskRegression_coxVariableName(object,
                                              model.frame = object.data)
    name.Xstrata <- coxInfo$stratavars.original
    ##terms.special <- prodlim::strip.terms(terms(coxFormula(object)), special = "strata")
    name.Xlp <- coxInfo$lpvars#attr(terms.special, "term.labels")
    name.X <- c(name.Xstrata,name.Xlp)
    only.strata <- if(length(name.Xlp)==0 && length(name.Xstrata)>0){TRUE}else{FALSE}
    if(only.strata){
        var.strata <- "strata"
        name.strata <- "strata"
    }else if(length(name.X)>0){
        var.strata <- "strata2"
        name.strata <- NULL
    }else{
        var.strata <- NULL
        name.strata <- NULL
    }
    
    ## 
    if(is.null(data)){
        originalData <- as.data.table(object.data)
        ## cannot use model.frame = TRUE otherwise the strata variable are combined into one
        ## add strata
        originalData[,c("strata") := riskRegression::coxStrata(object,
                                                               data = originalData,
                                                               sterms = coxInfo$strata.sterms,
                                                               strata.vars = coxInfo$stratavars,
                                                               strata.levels = coxInfo$strata.levels
                                                               )]
  
    }else{
        originalData <- copy(as.data.table(data))
    }

    ## initialize newdata
    if(is.null(newdata)){
        newdata <- originalData[,.SD,.SDcols=c("stop","strata",name.X)]
        init.newdata <- TRUE
    }else{        
        init.newdata <- FALSE
    }

    ## add prediction times and restrict to one observation per strata
    if(init.newdata){ 
        newdata[,c("Utimes") := .(.(sort(unique(originalData$stop)))), by = "strata"]
        newdata <- newdata[, .SD[1], by = name.X]

        if(identical(var.strata,"strata2")){
            newdata.X <- newdata[,.SD,.SDcols = name.X]
            tmp <- data.frame(lapply(1:NCOL(newdata.X),function(j){
                paste0(name.X[j],"=",newdata.X[[name.X[j]]])
            }))
            newdata[, c("strata2") := apply(tmp,1,paste0,collapse = ", ")]
        }
    }else{
        newdata[, c("Utimes") := .(.(sort(unique(originalData$stop))))]
        if(length(name.Xstrata)==0){
            newdata[, c("strata") := unique(originalData$strata)]
        }
        if(only.strata==FALSE){
            newdata[,c("strata2") := paste0("obs",1:.N)]
        }
    }

  ## compute predictions
  n.newdata <- NROW(newdata)    
  predC <- NULL    
  for(iObs in 1:n.newdata){ # iObs <- 1

      outRR <- riskRegression::predictCox(object,
                                          newdata = newdata[iObs],
                                          times = newdata[iObs,.SD$Utimes[[1]]],
                                          se = confint,
                                          keep.strata = FALSE,
                                          keep.newdata = TRUE,
                                          type = "survival")

      outRR$newdata <- newdata[iObs,.SD,.SDcols = union(c(name.X,"strata"),var.strata)]
      rrPred  <- as.data.table(outRR)
      
      ## add the number of events
      dt.tempo <- originalData[strata == newdata[iObs,strata],
                               .("n.censor" = sum(.SD$status == 0), "n.event" = sum(.SD$status == 1)),
                               by = "stop"]
      index.missing <- which(rrPred$times %in% dt.tempo$stop == FALSE)
      if(length(index.missing)>0){
          dt.tempo <- rbind(dt.tempo,
                            cbind(stop = rrPred$times[index.missing], n.censor = 0, n.event = 0))
      }

      dt.tempo <- merge(x = rrPred, y = dt.tempo, by.x = "times", by.y = "stop")
      predC <- rbind(predC, dt.tempo)
  }
    
    ## normalize for export
    if(confint==FALSE){
        predC[, c("survival.lower") := NA]
        predC[, c("survival.upper") := NA]
    }

    ## export
    res <- ggSurv(object = predC,
                  var.time = "times", var.survival = "survival",
                  confint = confint,
                  var.ci.inf = "survival.lower",
                  var.ci.sup = "survival.upper",
                  var.event = "n.event", var.censor = "n.censor", var.strata = var.strata, name.strata = name.strata, 
                  ...)
  
  return(invisible(res))
  
}
# }}}

## * ggSurv.data.table (code)
#' @rdname ggSurv
#' @export
ggSurv.data.table <- function(object, format = "data.table",
                              var.time = "time", var.survival = "survival", var.ci.inf = NULL, var.ci.sup = NULL, 
                              var.event = NULL, var.censor = NULL,  name.strata = NULL, var.strata = NULL, 
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
      gg.base <- ggplot2::ggplot(data = dt2.ggplot,
                                 mapping = aes(x = time, y = survival))
    }else{
      gg.base <- ggplot2::ggplot(data = dt2.ggplot,
                                 mapping = aes_string(x = "time", y = "survival",
                                                      color = var.strata, fill = var.strata))
    }
    
    # ci
    if(confint){
      if (is.null(var.strata)) {
        gg.base <- gg.base + ggplot2::geom_ribbon(aes(x = time, ymin = ci.inf, ymax = ci.sup, alpha = alpha.CIfill))
      }else{
        gg.base <- gg.base + ggplot2::geom_ribbon(data = dt2.ggplot,
                                                  aes(x = time, ymin = ci.inf, ymax = ci.sup),
                                                  alpha = alpha.CIfill)
        gg.base <- gg.base + ggplot2::geom_line(data = dt2.ggplot,
                                                aes(x = time, y = ci.inf),
                                                alpha = alpha.CIline)
        gg.base <- gg.base + ggplot2::geom_line(data = dt2.ggplot,
                                                aes(x = time, y = ci.sup),
                                                alpha = alpha.CIline)
      }
    }
    
                                        # step function
    gg.base <- gg.base + ggplot2::geom_line(data = dt2.ggplot,
                                            size = line.size)
    if (is.null(var.strata)) {
        gg.base <- gg.base + ggplot2::guides(fill = guide_legend(title = name.strata, title.position = legend.position),
                                             color = guide_legend(title = name.strata, title.position = legend.position)) #+ scale_colour_discrete(name = varStrata)
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
      gg.base <- gg.base + ggplot2::scale_shape_manual( values = values.ShapeLegend )
      gg.base <- gg.base + ggplot2::guides(shape = guide_legend(title = "Observation", title.position = legend.position ))
    }
    
    if(!is.null(title)){gg.base <- gg.base + ggplot2::ggtitle(title)}
    if(!is.null(text.size)){gg.base <- gg.base +  ggplot2::theme(text = element_text(size = text.size))}
    if(!is.null(ylim)){gg.base <- gg.base + ggplot2::coord_cartesian(ylim = ylim)}

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




