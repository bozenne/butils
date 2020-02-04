### descriptiveTable.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov  1 2018 (14:00) 
## Version: 
## Last-Updated: jan 31 2020 (11:50) 
##           By: Brice Ozenne
##     Update #: 286
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * descriptiveTable - doc
##' @title Descriptive Table for Continuous, Categorical, and Date Variables
##' @description Descriptive table for continuous, categorical, and date variables.
##' @name descriptiveTable
##' 
##' @param formula [formula] A formula where the left hand side (lhs) describe potential subgroups
##' and the right hand side (rhs) the variable to be described.
##' @param data [data.frame/data.table] Dataset
##' @param guess.categorical [integer,>0] When the type of the variables is not specified, numeric variables with
##' number of unique values lower than the specified value are treated as categorical variabels.
##' @param add.groupAll [logical] should the descriptives of the whole dataset be computed?
##' @param add.groupNA [logical] should the descriptives of the observations where the group variable is missing be computed?
##' @param add.groupVariable [logical] should the name of the variable used to generate the sub-groups be printed in the table.
##' @param FCT.center [function] function used to assess the center of the distribution.
##' @param FCT.spread [function] function used to assess the spread of the distribution.
##' @param test.categorical [function] function used to compute the p-value when comparing the categorical outcomes across sub-groups.
##' @param test.continuous [function] function used to compute the p-value when comparing the continuous outcomes across sub-groups.
##' 
##' @examples
##' data(veteran, package = "survival")
##'
##' descriptiveTable( ~ celltype + time + status +karno,
##'                  data = veteran)
##' descriptiveTable(trt ~ celltype + time + status +karno,
##'                  data = veteran)
##' tt <- descriptiveTable(trt ~ celltype + time + status +karno,
##'                        data = veteran, test.continuous = t.test)

## * descriptiveTable - code
##' @rdname descriptiveTable
##' @export
descriptiveTable <- function(formula, data, guess.categorical = 5,
                             add.groupAll = TRUE, add.groupNA = TRUE, add.groupVariable = FALSE,
                             FCT.center = "mean", FCT.spread = "sd",
                             test.categorical = NULL, test.continuous = NULL){

    data <- data.table::as.data.table(data)
    name.data <- names(data)

    ## ** name of the group variables
    formula.left <- update(formula,"~0")
    group.var <- all.vars(formula.left)

    if(length(group.var)> 0 && any(group.var %in% name.data == FALSE)){
        txt <- group.var[group.var %in% name.data == FALSE]
        stop("Argument \'data\' does not contain all variables in the lhs of the formula \n",
             "missing variable(s): \"",paste(txt, collapse = "\" \""),"\" \n")
    }
    if(length(group.var)> 0 && any(duplicated(group.var))){
        txt <- unique(group.var[duplicated(group.var)])
        stop("Incorrect specification of the formula, duplicated variable in the lhs \n",
             "missing variable(s): \"",paste(txt, collapse = "\" \""),"\" \n")
    }
    
    ## ** name of the descriptive variables
    formula.right <- update(formula,"0~.")
    formula.right.char <- as.character(formula.right)[3]
    X.var <- all.vars(formula.right)
    n.X <- length(X.var)

    if(any(X.var %in% name.data == FALSE)){
        txt <- X.var[X.var %in% name.data == FALSE]
        stop("Argument \'data\' does not contain all variables in the rhs of the formula \n",
             "missing variable(s): \"",paste(txt, collapse = "\" \""),"\" \n")
    }
    if(any(duplicated(X.var))){
        txt <- unique(X.var[duplicated(X.var)])
        stop("Incorrect specification of the formula, duplicated variable in the lhs \n",
             "missing variable(s): \"",paste(txt, collapse = "\" \""),"\" \n")
    }
    if(length(group.var)> 0 && length(intersect(X.var,group.var))>0){
        txt <- intersect(X.var,group.var)
        stop("Incorrect specification of the formula, one or more variables is both in the lhs and in the rhs of the formula \n",
             "variable(s): \"",paste(txt, collapse = "\" \""),"\" \n")
    }
    
    ## ** type of the descriptive variables
    type <- rep(NA, n.X)    
    vec.Xform <- strsplit(formula.right.char, split = "+", fixed = TRUE)[[1]]
    vec.Xform <- trimws(vec.Xform, which = "both") 

    index.type <- grep("(",vec.Xform, fixed = TRUE)
    for(iX in 1:n.X){
        if(iX %in% index.type){            
            type[iX] <- strsplit(vec.Xform[iX], split = "(", fixed = TRUE)[[1]][1]
        }else{
            type[iX] <- .guessType(data[[X.var[iX]]], guess.categorical = guess.categorical)
        }
    }
    
    type <- trimws(type, which = "both") ## remove white spaces

    ## ** define groups of observations
    group <- "All"
    nObs.group <- NROW(data)
    if(length(group.var)>0){
        dt.tempo <- data[,.N, by = group.var]
        if(add.groupVariable){
            name.subgroups <- paste(paste(group.var,collapse="/"),"=",dt.tempo[, apply(.SD,1,paste,collapse="/") ,.SDcols = group.var])
        }else{
            name.subgroups <- paste(dt.tempo[, apply(.SD,1,paste,collapse="/") ,.SDcols = group.var])
        }
        group <- c(group, name.subgroups)
        nObs.group <- c(nObs.group,
                        dt.tempo[["N"]])
    }
    name.group <- paste0(group, " (n=",nObs.group,")")

    ## ** all observations    
    ls.res <- lapply(1:n.X, function(iX){
        iOut <- .summaryStat(values = data[[X.var[iX]]], type = type[iX],
                             FCT.center = FCT.center, FCT.spread = FCT.spread)
        return(cbind(group = name.group[1], variable = X.var[iX], iOut))
    })

    ## ** sub-groups
    if(length(group.var)>0){
        iData <- data[,list(data = list(.SD)), by = group.var]
        n.group <- NROW(iData)
        
        for(iX in 1:n.X){
            for(iG in 1:n.group){
                iOut <- .summaryStat(values = iData[[2]][[iG]][[X.var[iX]]], type = type[iX],
                                     FCT.center = FCT.center, FCT.spread = FCT.spread)
                ls.res[[iX]] <- rbind(ls.res[[iX]],cbind(group = name.group[iG+1], variable = X.var[iX], iOut))

            }
        }

        if(!is.null(test.categorical) && any(c("categorical","binary") %in% type)){
            index.categorical <- which(type %in% c("categorical","binary"))
            n.categorical <- length(index.categorical)
            p.value.categorical <- vector(mode = "list", length = n.categorical)
            names(p.value.categorical) <- X.var[index.categorical]
                
            for(iX in 1:n.categorical) { ## iX <- 1
                lsXY <- lapply(iData[[2]],"[[",X.var[index.categorical[iX]])

                if(any(sapply(lsXY,length)==0) || any(sapply(lsXY, function(iXX){sum(!is.na(iXX))})==0)){
                    p.value.categorical[[iX]] <- NA
                }else{
                    n.group <- length(lsXY)
                    df.testCategorical <- data.frame(x = unlist(lapply(1:n.group, function(iG){rep(iG, length(lsXY[[iG]]))})),
                                                     y = unlist(lsXY))                    

                    p.value.categorical[[iX]] <- try(do.call(test.categorical, args = as.list(df.testCategorical))$p.value, silent = FALSE)
                }
            }
        }else{
            test.categorical <- NULL
        }

        if(!is.null(test.continuous) && any("continuous" %in% type)){
            index.continuous <- which(type %in% c("continuous"))
            n.continuous <- length(index.continuous)
            p.value.continuous <- vector(mode = "list", length = n.continuous)
            names(p.value.continuous) <- X.var[index.continuous]

            for(iX in 1:n.continuous) { ## iX <- 1
                lsXY <- lapply(iData[[2]],"[[",X.var[index.continuous[iX]])
                
                if(any(sapply(lsXY,length)==0) || any(sapply(lsXY, function(iXX){sum(!is.na(iXX))})==0)){
                    p.value.continuous[[iX]] <- NA
                }else{
                    n.group <- length(lsXY)
                    df.testContinuous <- data.frame(group = unlist(lapply(1:n.group, function(iG){rep(iG, length(lsXY[[iG]]))})),
                                                    value = unlist(lsXY))
                    p.value.continuous[[iX]] <- try(do.call(test.continuous, args = list(formula = value ~ group, data = df.testContinuous))$p.value, silent = FALSE)
                }
            }
        }else{
            test.continuous <- NULL
        }
    }else{
        test.categorical <- NULL
        test.continuous <- NULL
    }
        
    ## ** merge
    type.merge <- type
    type.merge[type.merge == "binary"] <- "categorical"
    out <- tapply(ls.res, type.merge, function(iList){
        iOut <- do.call(rbind, iList)
        
        name.rm <- NULL
        if(add.groupAll == FALSE){ ## iList <- ls.res[[1]]
            name.rm <- c(name.rm,name.group[which(group=="All")])
        }
        if(add.groupNA == FALSE){ ## iList <- ls.res[[1]]
            name.rm <- c(name.rm,name.group[which(group=="NA")])
        }
        if(length(name.rm)>0){
            iOut <- iOut[iOut$group != name.rm,]
        }
        return(iOut)

    })

    if(!is.null(test.categorical)){
        attr(out$categorical,"p.value") <- p.value.categorical
    }
    if(!is.null(test.continuous)){
        attr(out$continuous,"p.value") <- p.value.continuous
    }
    
    ## ** export
    class(out) <- "descriptiveTable"
    attr(out, "name.group") <- name.group
    attr(out, "FCT.center") <- FCT.center
    attr(out, "FCT.spread") <- FCT.spread
    return(out)

}

## * print.descriptiveTable
##' @title Print Function for the Descriptive Table
##' @description Print function for the descriptive table.
##'
##' @param x output of \code{descriptiveTable}.
##' @param print [logical] should the descriptive table be printed?
##' @param print.minmax [logical] should minimum and maximum values be displayed?
##' @param print.spread [logical] should the spread of the values be displayed?
##' @param digit.frequency [integer, >=0] number of digit when printing frequencies.
##' @param digit.center [integer, >=0] number of digit when printing center parameters.
##' @param digit.spread [integer, >=0] number of digit when printing spread parameters.
##' @param digit.p [integer, >=0] number of digit when printing p-values.
##' @param format.date [character] the format used to output dates.
##' @param ... Not used. For compatibility with the generic method.
##' 
##' @export
print.descriptiveTable <- function(x, print = TRUE,
                                   print.minmax = TRUE,
                                   print.spread = TRUE,
                                   digit.frequency = 2,
                                   digit.center = 2,
                                   digit.spread = 2,
                                   digit.p = 2,
                                   format.date = "%Y-%m-%d",
                                   ...){

    name.group <- attr(x, "name.group")
    n.group <- length(name.group)
    name.type <- names(x)
    n.type <- length(name.type)

    out <- vector(mode = "list", length = n.type)
    names(out) <- name.type

    ## ** rounding
    if("date" %in% name.type){
        vec.value <- c("center","[min;max]","n.NA")
        if(print.minmax == FALSE){
            vec.value <- vec.value[vec.value!="[min;max]"]
        }

        out$date <- data.table::copy(data.table::as.data.table(x[["date"]]))
        out$date[, c("center") := format(.SD[["center"]], format = format.date)]
        out$date[, c("min") := format(.SD[["min"]], format = format.date)]
        out$date[, c("max") := format(.SD[["max"]], format = format.date)]
        out$date[,c("[min;max]") := paste(" [",.SD$min,";",.SD$max,"]",sep="")]
        out$date <- dcast(out$date, value.var = vec.value, formula = variable ~ group)
    }
    if("constant" %in% name.type){
        out$constant <- data.table::copy(data.table::as.data.table(x[["constant"]]))
        out$constant[, c("level") := as.character(.SD$level)]
        out$constant <- dcast(out$constant, value.var = c("n","frequency","n.NA"), formula = variable + level ~ group)        
    }
    if("categorical" %in% name.type){
        out$categorical <- data.table::copy(data.table::as.data.table(x[["categorical"]]))
        out$categorical[, c("frequency") := round(100*.SD[["frequency"]], digits = digit.frequency)]
        out$categorical[, c("level") := as.character(.SD$level)]
        out$categorical <- dcast(out$categorical, value.var = c("n","frequency"), formula = variable + level ~ group)

        p.value <- attr(x[["categorical"]],"p.value")
        if(!is.null(p.value)){
            p.value <- unlist(lapply(p.value, function(iP){
                if(is.numeric(iP)){
                    return(format.pval(iP, digits=digit.p,eps=10^{-digit.p}))
                }else{
                    return(NA)
                }
            }))
            out$categorical[, p.value :=  c(rep(NA, times=.N-1), p.value[.GRP]), by="variable"]
        }
    }
    if("continuous" %in% name.type){
        vec.value <- c("center","spread","[min;max]","n.NA")
        if(print.spread == FALSE){
            vec.value <- vec.value[vec.value!="spread"]
        }
        if(print.minmax == FALSE){
            vec.value <- vec.value[vec.value!="[min;max]"]
        }
               
        out$continuous <- data.table::copy(data.table::as.data.table(x[["continuous"]]))
        out$continuous[, c("center") := round(.SD[["center"]], digits = digit.center)]
        out$continuous[, c("min") := round(.SD[["min"]], digits = digit.center)]
        out$continuous[, c("max") := round(.SD[["max"]], digits = digit.center)]
        out$continuous[, c("spread") := paste0("(",round(.SD[["spread"]], digits = digit.center),")")]
        out$continuous[, c("[min;max]") := paste0("[",.SD[["min"]],";",.SD[["max"]],"]")]
        out$continuous <- dcast(out$continuous, value.var = vec.value, formula = variable ~ group)

        p.value <- attr(x[["continuous"]],"p.value")
        if(!is.null(p.value)){
            p.value <- unlist(lapply(p.value, function(iP){
                if(is.numeric(iP)){
                    return(format.pval(iP, digits=digit.p,eps=10^{-digit.p}))
                }else{
                    return(NA)
                }
            }))
            out$continuous[, c("p.value") :=  p.value]
        }
    }

    ## ** rename according to center and scale
    if(attr(x,"FCT.center") %in% c("mean","median")){
        if("date" %in% name.type){
            new.names <- gsub("center",attr(x,"FCT.center"),names(out$date), fixed = TRUE)
            setnames(out$date, old = names(out$date), new = new.names)
        }
        if("continuous" %in% name.type){
            new.names <- gsub("center",attr(x,"FCT.center"),names(out$continuous), fixed = TRUE)
            setnames(out$continuous, old = names(out$continuous), new = new.names)
        }
    }
    if(attr(x,"FCT.spread") %in% c("sd","IQR")){
        if("continuous" %in% name.type){
            new.names <- gsub("spread",attr(x,"FCT.spread"),names(out$continuous), fixed = TRUE)
            setnames(out$continuous, old = names(out$continuous), new = new.names)
        }
    }

    ## ** print
    out.print <- vector(mode = "list", length = n.type)
    if(print){

        for(iType in 1:n.type){
            iNames.out <- names(out[[iType]])

            iOut.print <- data.table::copy(out[[iType]])
            iOut.print[,c("variable") := as.character(.SD$variable)]
            iOut.print[,c("variable") := c(.SD$variable[1],rep("",.N-1)), by = "variable", .SDcols = "variable"]

            ## re-order columns
            if(name.type[iType] %in% c("binary","categorical")){
                iNew.order <- c("variable","level",iNames.out[unlist(lapply(name.group, grep, iNames.out, fixed = TRUE))])
            }else{
                iNew.order <- c("variable",iNames.out[unlist(lapply(name.group, grep, iNames.out, fixed = TRUE))])
            }
            setcolorder(iOut.print, neworder = iNew.order)

            ## rename columns
            toRemove <- paste0("_",name.group)
            iColnames <- names(iOut.print)
            for(iRm in toRemove){
                iColnames <- gsub(iRm,"",iColnames,fixed = TRUE)
            }

            iOut.print2 <- rbind(iColnames,
                                 as.data.frame(iOut.print))

            ## rename first header
            names(iOut.print2)[names(iOut.print2) %in% c("variable","level","p.value")] <- ""
            for(iG in name.group){
                names(iOut.print2)[grep(iG, names(iOut.print2), fixed = TRUE)] <- paste0(" ",iG)
            }
            names(iOut.print2)[duplicated(names(iOut.print2))] <- ""

            cat("    ** ",name.type[iType]," variables ** \n", sep = "")
            print(iOut.print2, row.names = FALSE)
            cat("\n")

            out.print[[iType]] <- iOut.print2
        }
    }

    officer.print <- lapply(out.print, function(iT){
        iT2 <- rbind(names(iT),iT)
        names(iT2) <- NULL
        iT2[is.na(iT2)] <- ""
        return(iT2)
    })
    
    return(invisible(list(table = out,
                          table.print = out.print,
                          table.officer = officer.print)))
}
    


## * .guessType
.guessType <- function(x, guess.categorical){

    n.unique <- length(unique(na.omit(x)))
    if(inherits(x, "Date") || inherits(x, "POSIXct") || inherits(x, "POSIXt")){
        return("date")
    }else if(n.unique == 1){
        return("constant")
    }else if(n.unique == 2){
        return("binary")
    }else if(n.unique <= guess.categorical || is.character(x) || is.factor(x)){
        return("categorical")
    }else if((is.numeric(x) || is.integer(x)) ){
        return("continuous")
    }else {
        stop("Did not managed to guess the type. Please specify the type in the formula \n")
    }
}

## * .summaryStat
.summaryStat <- function(values, type, FCT.center, FCT.spread){

    if(type == "date"){
        if(all(is.na(values))){
            out <- data.frame(n = length(values),
                              n.NA = sum(is.na(values)),                 
                              center = NA,
                              spread = NA,
                              min = NA,
                              max = NA,
                              stringsAsFactors = FALSE)
        }else{
            out <- data.frame(n = length(values),
                              n.NA = sum(is.na(values)),                 
                              center = do.call(FCT.center, args = list(values, na.rm = TRUE)),
                              spread = do.call(FCT.spread, args = list(values, na.rm = TRUE))/(3600*24),
                              min = min(values, na.rm = TRUE),
                              max = max(values, na.rm = TRUE),
                              stringsAsFactors = FALSE)
        }
    }else if(type == "constant"){
        out <- data.frame(n = length(values),
                          n.NA = sum(is.na(values)),
                          level = values[1],
                          frequency = 1,
                          stringsAsFactors = FALSE)
    }else if(type %in% c("binary","categorical")){
        table.values <- table(values, useNA = "ifany")
        out <- data.frame(n = as.double(table.values),
                          level = names(table.values),
                          frequency = as.double(table.values)/sum(table.values),
                          stringsAsFactors = FALSE)
    }else if(type == "continuous"){
        if(all(is.na(values))){
            out <- data.frame(n = length(values),
                              n.NA = sum(is.na(values)),
                              center = NA,
                              spread = NA,
                              min = NA,
                              max = NA,
                              stringsAsFactors = FALSE)
        }else{
            out <- data.frame(n = length(values),
                              n.NA = sum(is.na(values)),
                              center = do.call(FCT.center, args = list(values, na.rm = TRUE)),
                              spread = do.call(FCT.spread, args = list(values, na.rm = TRUE)),
                              min = min(values, na.rm = TRUE),
                              max = max(values, na.rm = TRUE),
                              stringsAsFactors = FALSE)
        }
    } 
    return(out)
}


######################################################################
### descriptiveTable.R ends here
