### print.detectOutlier.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt  5 2017 (10:15) 
## Version: 
## last-updated: okt  5 2017 (10:15) 
##           By: Brice Ozenne
##     Update #: 4
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * print.detectOutlier
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

## * print.ls_detectOutlier
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



#----------------------------------------------------------------------
### print.detectOutlier.R ends here
