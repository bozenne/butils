# {{{ scanOutlier
#' @title Search potential outliers
#' @description Search potential outliers in a dataset
#' 
#' @param data the dataset
#' @param id the name of the column containg a unique identifier for each observation
#' @param method.numeric the method used to detect outliers in a numeric variable
#' @param method.factor the method used to detect outliers in a categorical variable
#' @param method.id the method used to detect outliers in identifiers
#' @param args.num arguments to be passed to the method.numeric 
#' @param args.factor arguments to be passed to the method.factor 
#' 
#' @export
scanOutlier <- function(data, id, 
                        method.numeric ="numOutlier", method.factor ="factorOutlier", method.id ="checkId",
                        args.num = NULL, args.factor = NULL){
  
  if(!data.table::is.data.table(data)){
    data <- data.table::as.data.table(data)
  }
  
  names.data <- names(data)
  test.missingID <- missing(id)
  
  output <- lapply(names.data, function(name){
    x <- data[[name]]
    
    if(!test.missingID && name == id){
      
      return(do.call(checkDuplicated, args = list(x)))
      
    }else if(is.numeric(x)){
      res <- checkUnique(x, test = TRUE)
      if(class(res)[1] == "detectOutlier"){
        return(res)
      }else{
        return(do.call(method.numeric, args = c(list(x), args.num)))    
      }
      
    }else if(is.factor(x)){
      
      res <- checkUnique(x, test = TRUE)
      if(class(res)[1] == "detectOutlier"){
        return(res)
      }else{
        return(do.call(method.factor, args = c(list(x), args.factor)))    
      }
      
    }else if(is.character(x)){
      
      return(checkUnique(x, test = FALSE))
      
    }else{
      
      return("unknow type")
      
    }
  })
  names(output) <- names.data
  
  class(output) <- "ls_detectOutlier"
  return(output)
}
# }}}

# {{{ check functions

# {{{ doc
#' @title Identify outlier
#' @name identifyOutlier
#' @description Identify numeric or factor outliers
#'
#' @param x a vector of numeric
#' @param na.rm should na be removed.
#' @param threshold the threshold for defining an outlier regarding the prevalence (categorical variables)
#' @param th.gaussian the threshold for defining an outlier when using the mean
#' @param th.hampel the threshold for defining an outlier when mad
#' @param th.boxplot the threshold for defining the IQR
#' @param type the type of robust metric for assessing the "average" value
#' @param use.NA argument passed to \code{table}
#' 
#' @references idea from http://www.r-bloggers.com/finding-outliers-in-numerical-data/
#'
#' @examples 
#' \dontrun{
#' numOutlier(rnorm(1e3))
#' }
# }}}

# {{{ numOutlier
#' @rdname identifyOutlier
numOutlier <- function(x, type = "auto",
                       th.gaussian = 3, th.hampel = 3, th.boxplot = 1.5, 
                       na.rm = FALSE){
  
  validCharacter(type, valid.values = c("gaussian", "hampel", "boxplot", "auto"), valid.length = 1)
  if(type == "auto"){
    if(stats::mad(x, na.rm = na.rm) == 0){type <- "boxplot"}else{type <- "hampel"}
  }
  
  if(any(is.infinite(x))){
    index.infinite <- which(is.infinite(x))
    output <- list(index = index.infinite,
                   value = x[index.infinite], 
                   type = "infinite",
                   details = list(),
                   discrepancy = rep("Inf", length(index.infinite)),
                   display = list())
    return(output)
  }
  
  threshold <- switch(type,
                      "gaussian" = th.gaussian,
                      "hampel" = th.hampel,
                      "boxplot" = th.boxplot)
  
  if(type == "gaussian"){
    center <- rep(mean(x, na.rm = na.rm), 2)
    scale <- sd(x, na.rm = na.rm)
  }else if(type == "hampel"){
    center <- rep(median(x, na.rm = na.rm), 2)
    scale <- stats::mad(x, na.rm = na.rm)
  }else if(type == "boxplot"){
    center <- quantile(x, probs = c(0.25,0.75), na.rm = na.rm)
    scale <- stats::IQR(x, na.rm = na.rm)
  }
  
  limitInf <- center[1] - threshold*scale
  limitSup <- center[2] + threshold*scale
  if(is.na(limitInf) || is.na(limitSup)){
    stop("detectOutlier: \"x\" contains NA values or insufficient number of values \n")
  }
  
  outliers <- union(which(x < limitInf),
                    which(x > limitSup))
  outliers.th <- apply(cbind((center[1]-x[outliers])/scale,-(center[2]-x[outliers])/scale), 1, 
                       function(x){x[x>0]})
  
  
  #### export
  output <- list(index = outliers,
                 type = "numeric",
                 value = x[outliers], 
                 details = list(type = type,
                                interval = c(limitInf, limitSup)),
                 discrepancy = outliers.th,
                 display = list(args = boxplot(x, plot = FALSE),
                                method = "boxplot"))
  
  class(output) <- "detectOutlier"
  return(output)
}
# }}}

# {{{ factorOutlier
#' @rdname identifyOutlier
factorOutlier <- function(x, threshold = 0.01, use.NA = "ifany"){
  
  groups <- unique(x)
  n.group <- length(groups)
  
  tabx <- table(x, useNA = use.NA)
  prevalence <- tabx/sum(tabx)
  test.outlier <- prevalence < threshold/n.group
  group.outlier <- names(tabx)[test.outlier]
  
  index.outlier <- which(x %in% group.outlier)
  output <- list(index = index.outlier,
                 value = x[index.outlier],
                 type = "factor",
                 details = list(threshold = threshold/n.group),
                 discrepancy = setNames(prevalence[test.outlier], group.outlier),
                 display = list(args = tabx,
                                method = "barplot"))
  
  class(output) <- "detectOutlier"
  return(output)
}
# }}}

# {{{ checkDuplicated
checkDuplicated <- function(x){
  
  index.duplicated <- which(duplicated(x))
  
  if(length(index.duplicated)>0){
    discrepancy <- unique(x[index.duplicated])
  }else{
    discrepancy <- integer(0L)
  }
  
  output <- list(index = index.duplicated,
                 value = x[index.duplicated],
                 type = "duplicated",
                 details = list(levels = unique(x)),
                 discrepancy = discrepancy,
                 display = list())
  
  class(output) <- "detectOutlier"
  return(output)
}
# }}}

# {{{ checkUnique
checkUnique <- function(x, test = FALSE){
  
  if(length(unique(x))==1){
    output <- list(index = 1:length(x),
                   value = x,
                   type = "unique",
                   details = list(levels = x[1]),
                   discrepancy = NA,
                   display = list())
    class(output) <- "detectOutlier"
    return(output)
    
  }else{
    
    if(test){
      return(FALSE)
    }else{
      output <- list(index = integer(0L),
                     value = x[integer(0L)],
                     type = "unique",
                     details = list(),
                     discrepancy = integer(0L),
                     display = list())  
      class(output) <- "detectOutlier"
      return(output)
    }
    
  }
  
}
# }}}

# }}}

# {{{ print 

# {{{ print.detectOutlier
print.detectOutlier <- function(x, type = "value", ...){
  
  validCharacter(type, valid.length = 1, valid.values = c("value","index"), method = "print.detectOutlier")
  
  if(x$type == "numeric"){
    cat.message <- "potential outliers: \n"
    vec <- setNames(x[[type]], x$discrepancy)
  }else if(x$type == "factor"){
    cat.message <- "levels with few instances: \n"
    vec <- x$discrepancy
  }else if(x$type == "duplicated"){
    cat.message <- "duplicated levels: \n"
    vec <- x$discrepancy
  }else if(x$type == "infinite"){
    cat.message <- NULL
    vec <- paste0("number of infinite values: ",length(x$discrepancy))
  }else if(x$type == "unique"){
    cat.message <- NULL
    vec <- paste0("only one level: ", x$details$levels)
  }
  if(length(x$index)==0){
    cat("no outlier \n")
  }else{
    if(!is.null(cat.message)){cat(cat.message)}
    print(vec)
  }
  
}
# }}}

# {{{ print.ls_detectOutlier
print.ls_detectOutlier <- function(x, only.outlier = TRUE, ...){
  
  names.x <- names(x)
  
  res <- lapply(names.x, function(name){
    if(length(x[[name]]$index)>0 || only.outlier == FALSE){
      cat("# ",name," (",x[[name]]$type,") - ", sep = "")
      print(x[[name]], ...)
      return(TRUE)
    }else{
      return(FALSE)
    }
  })
  return(invisible())
}
# }}}

# }}}

# {{{ plot 

# {{{ plot.detectOutlier
#' @export
plot.detectOutlier <- function(x, ...){
  
  if(length(x$display)>0){
    do.call(x$display$method, args = x$display$args)
  }else{
    message("no available plot \n")
  }
  
}
# }}}

# }}}
