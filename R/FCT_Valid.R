#**********************************************************************
#**********************************************************************
#*************         validation functions         *******************
#**********************************************************************
#**********************************************************************
#

#' @name validFCTs
#' @aliases validClass
#' @aliases validDimension
#' @aliases validInteger
#' @aliases validLogical
#' @aliases validNames
#' @aliases validNumeric
#' @aliases validPath
#' 
#' @title check the validity of the arguments given by the user
#' 
#' @param value1
#' @param value2
#' @param name1
#' @param name2
#' @param validClass
#' @param validDimension
#' @param validLength
#' @param validValues
#' @param refuse.NULL
#' @param refuse.duplicates
#' @param superClasses
#' @param type
#' @param method
#' @param addPP
#' 
#' @rdname validFCTs
#' @export
validCharacter <- function(value1, name1 = as.character(substitute(value1)), validLength, 
                           validValues = "character", refuse.NULL = TRUE, refuse.duplicates = FALSE, 
                           method = NULL, addPP = TRUE){
  
  if(!is.null(method) && addPP){
    method <- paste0(method, ": ")
  }
  
  if(is.null(value1)){
    
    if(refuse.NULL == TRUE){
      stop(method, "wrong specification of \'", name1, "\' \n", 
           "\'", name1, "\' must not be NULL \n")
    }
    
  }else{
    
    #### check size
    n.value1 <- length(value1)
    
    if(!is.null(validLength) && n.value1 %in% validLength == FALSE){
      stop(method, "wrong specification of \'", name1, "\' \n", 
           "\'", name1, "\' must have length ", paste(validLength, collapse = " or "), "  \n", 
           "length(", name1, ") : ", n.value1, "\n")
    }
    
    #### check duplicates
    if(refuse.duplicates == TRUE && any(duplicated(value1))){
      stop(method, "wrong specification of \'", name1, "\' \n", 
           "\'", name1, "\' contains duplicated values : \n",        
           "\"",paste(unique(value1[duplicated(value1)]), collapse = "\" \""), "\" \n")
    }
    
    #### check values
    if(identical(validValues,"character")){
      
      if(any(is.character(value1) == FALSE)){
        stop(method, "wrong specification of \'", name1, "\' \n", 
             "\'", name1, "\' must be a ", if(n.value1 == 1){"character"}else{"vector of characters"}," \n", 
             "is(", name1, ") : ", paste(is(value1), collapse = " "), "\n")
      }
      
    } else if(identical(validValues,"character_or_logical")){
      
      if(any( (is.character(value1) == FALSE) * (is.logical(value1) == FALSE) > 0 )){
        stop(method, "wrong specification of \'", name1, "\' \n", 
             "\'", name1, "\' must be a ", if(n.value1 == 1){"character or logical"}else{"vector of characters or logicals"}," \n", 
             "is(", name1, ") : ", paste(is(value1), collapse = " "), "\n")
      }
      
    } else if(!is.null(validValues) && any(value1 %in% validValues == FALSE)){
      
      stop(method, "wrong specification of \'", name1, "\' \n", 
           "valid values for \'", name1, "\' : ", if(refuse.NULL == FALSE){"NULL"}, " \"", paste(validValues, collapse = "\" \""), "\" \n", 
           "refused value",if(sum(value1 %in% validValues == FALSE)>1){"s"}," for \'", name1, "\' : \"", paste(value1[value1 %in% validValues == FALSE], collapse = "\" \""), "\"\n")
      
    }
    
  }
  
  return(invisible(TRUE))
  
}

#' @rdname validFCTs
#' @export
validClass <- function(value1, name1 = as.character(substitute(value1)), validClass, 
                       superClasses = TRUE, method = NULL, addPP = TRUE){
  
  if(!is.null(method) && addPP){
    method <- paste0(method, ": ")
  }
  
  if(superClasses == TRUE){
    
    if( all(is(value1) %in% validClass == FALSE) ){
      stop(method, "wrong specification of \'", name1, "\' \n", 
           "superclass of \'", name1, "\' must be one of the following \"", paste(validClass,collapse="\" \""), "\"  \n", 
           "proposed superclass : \"", paste(is(value1),collapse="\" \""), "\" \n")
    }  
    
  }else{
 
    if( class(value1) %in% validClass == FALSE){
      stop(method, "wrong specification of \'", name1, "\' \n", 
           "class of \'", name1, "\' must be \"", paste(validClass,collapse="\" \""),"\"  \n", 
           "proposed class : ", class(value1)[[1]], "\n")
    }  
    
  }
  
  return(invisible(TRUE))
  
}

#' @rdname validFCTs
#' @export
validDimension <- function(value1, value2 = NULL, name1 = as.character(substitute(value1)), name2 = as.character(substitute(value2)),
                           validDimension = NULL,
                           type = c("NROW","NCOL"), method = NULL, addPP = TRUE){
  
  if(!is.null(method) && addPP){
    method <- paste0(method, ": ")
  }
  
  n.type <- length(type)
  
  #### dimension 1
  testDimension <- sapply(1:n.type, function(x){
    do.call(type[x], list(value1))
  })
  
  #### dimension 2
  
  
  if(is.null(validDimension)){
    
    validDimension <- sapply(1:n.type, function(x){
      do.call(type[x], list(value2))
    })
    test.validDimension <- TRUE
    
  }else if(is.null(name2)){
    
    test.validDimension <- FALSE
    
  }else{
    
    test.validDimension <- TRUE
    
  }
  
  #### main
  for(iter_type in 1:n.type){
    
    if(testDimension[iter_type] != validDimension[iter_type]){
      
      if(test.validDimension){
        stop(method, "dimension mismatch between argument \'", name1, "\' and argument \'", name2, "\' \n", 
             type[iter_type],"(", name1, ") = ", testDimension[iter_type], " \n", 
             type[iter_type],"(", name2, ") = ", validDimension[iter_type], " \n")  
      }else{
        stop(method, "dimension mismatch between argument \'", name1, "\' and argument \'", name2, "\' \n", 
             type[iter_type],"(", name1, ") = ", testDimension[iter_type], " \n", 
             type[iter_type],"(", name2, ") = ", validDimension[iter_type], " \n")
        
      }
      
    }
    
  }
    
  return(invisible(TRUE))
}
  
#' @rdname validFCTs
#' @export
validInteger <- function(value1, name1 = as.character(substitute(value1)), validLength, 
                         validValues = NULL, min = NULL, max = NULL, 
                         refuse.NA = TRUE, refuse.NULL = TRUE, refuse.duplicates = FALSE, 
                         method = NULL, addPP = TRUE){
  
  if(!is.null(method) && addPP){
    method <- paste0(method, ": ")
  }
  
  validNumeric(value1 = value1, name1 = name1, validLength = validLength, min = min, max = max, 
               refuse.NA = refuse.NA, refuse.NULL = refuse.NULL, refuse.duplicates = refuse.duplicates, method = method)
  
  #### check integer
  if(!is.null(value1) && any(value1 %% 1 > 0)){
    stop(method, "wrong specification of \'", name1, "\' \n", 
         "\'", name1, "\' must contain integers not doubles \n",        
         "invalid value(s) in ", name1, " : ", paste(value1[value1 %% 1 > 0], collapse = " "), "\n")
  }
  
  return(invisible(TRUE))
}

#' @rdname validFCTs
#' @export
validLogical <- function(value1, name1 = as.character(substitute(value1)), validLength, 
                         refuse.NULL = TRUE, refuse.NA = TRUE, 
                         method = NULL, addPP = TRUE){
  
  if(!is.null(method) && addPP){
    method <- paste0(method, ": ")
  }
  
  if(is.null(value1)){
    
    #### NULL
    if(refuse.NULL == TRUE){
      stop(method, "wrong specification of \'", name1, "\' \n", 
           "\'", name1, "\' must be logical ",if(refuse.NA == FALSE){"or NA"}," and not NULL \n")
    }
    
  }else{ 
    
    #### Size
    if(!is.null(validLength) && length(value1) %in% validLength == FALSE){
      stop(method, "wrong specification of \'", name1, "\' \n", 
           "\'", name1, "\' must have length ", paste(validLength, collapse = " or "), "  \n", 
           "length(", name1, ") : ", length(value1), "\n")
    } 
    
    #### Type
    if(any(is.logical(value1) == FALSE)){
      stop(method, "wrong specification of \'", name1, "\' \n", 
           "\'", name1, "\' must be ", if(refuse.NULL == FALSE){"NULL or "}, if(refuse.NA == FALSE){"NA or "},"TRUE or FALSE \n",        
           "is(", name1, ") : ", paste(is(value1), collapse = " "), "\n")
    }
    
    if(refuse.NA == TRUE && any(is.na(value1)) ){
      stop(method, "wrong specification of \'", name1, "\' \n", 
           "\'", name1, "\' must be logical ",if(refuse.NULL == FALSE){"or NULL"}," and not NA \n")
    }
    
  }
  
  return(invisible(TRUE))
}

#' @rdname validFCTs
#' @export
validNames <- function(value1, name1 = as.character(substitute(value1)), refuse.NULL = TRUE,
                       validLength = NULL, validValues = NULL, requiredValues = NULL, forbiddenValues = NULL,
                       method = NULL, addPP = TRUE){
  
  if(!is.null(method) && addPP){
    method <- paste0(method, ": ")
  }
  
  ## type
  if(is.matrix(value1)){
    value1 <- colnames(value1)
  }
  
  if(is.data.table(value1) || is.data.frame(value1) || is.list(value1)){
    value1 <- names(value1)
  }
  
  ## tests
  if(is.null(value1)){
    
    if(refuse.NULL == TRUE){
    stop(method, "wrong specification of \'", name1, "\' \n", 
         "names of \'", name1, "\' must not be NULL \n")
    }
    
  }else{
    
    #### check size
    n.value1 <- length(value1)
    
    if(!is.null(validLength) && n.value1 %in% validLength == FALSE){
      stop(method, "wrong specification of \'", name1, "\' \n", 
           "\'", name1, "\' must have ", paste(validLength, collapse = " or ")," names  \n", 
           "length(names(", name1, ")) : ", n.value1, "\n")
    }
    
    #### check content
    
    if(!is.null(requiredValues) && any(requiredValues %in% value1 == FALSE)){
      
      stop(method, "wrong specification of \'", name1, "\' \n", 
           "\'", name1, "\' must contains specific names \n",
           "missing names : \"",paste(requiredValues[requiredValues %in% value1 == FALSE], collapse = "\" \""),"\" \n", 
           "proposed names : \"", paste(value1, collapse = "\" \""), "\"\n")  
      
    }
    
    if(!is.null(validValues) && any(value1 %in% validValues == FALSE)){
      
      stop(method, "wrong specification of \'", name1, "\' \n", 
           "valid names for \'", name1, "\' : \"",paste(validValues, collapse = "\" \""),"\" \n", 
           "refused names : \"", paste(value1[value1 %in% validValues == FALSE], collapse = " "), "\"\n")  
      
    }
    
    if(!is.null(forbiddenValues) && any(value1 %in% forbiddenValues)){
      
      stop(method, "wrong specification of \'", name1, "\' \n", 
           "forbidden names for \'", name1, "\' : \"",paste(forbiddenValues, collapse = "\" \""),"\" \n", 
           "refused names : \"", paste(value1[value1 %in% forbiddenValues], collapse = " "), "\"\n")  
      
    }
    
    if(any(duplicated(value1))){
      stop(method, "wrong specification of \'", name1, "\' \n", 
           name1, " must not contain duplicated names \n", 
           "duplicated names : \"", paste(value1[duplicated(value1)], collapse = " "), "\"\n")  
    }
    
  }
  
  return(invisible(TRUE))
  
}
#' @rdname validFCTs
#' @export
validNumeric <- function(value1, name1 = as.character(substitute(value1)), validLength,
                         validValues = NULL , min = NULL, max = NULL,
                         refuse.NA = TRUE, refuse.NULL = TRUE, refuse.duplicates = FALSE, 
                         method = NULL, addPP = TRUE){
  
  if(!is.null(method) && addPP){
    method <- paste0(method, ": ")
  }
  
  if(is.null(value1)){
    
    if(refuse.NULL == TRUE){
      stop(method, "wrong specification of \'", name1, "\' \n", 
           "\'", name1, "\' must not be NULL \n")
    }
    
  }else{
    
    #### check length
    if(!is.null(validLength) && length(value1) %in% validLength == FALSE){
      stop(method, "wrong specification of \'", name1, "\' \n", 
           "\'", name1, "\' must have length ", paste(validLength, collapse = " or "), "  \n", 
           "length(", name1, ") : ", length(value1), "\n")
    }
    
    #### check NA
    if(refuse.NA == TRUE && any(is.na(value1))){
      stop(method, "wrong specification of \'", name1, "\' \n", 
           "\'", name1, "\' must not contain NA \n", 
           "index of NA values : ", which(paste(is.na(value1), collapse = " ")), "\n")
    }
    
    #### check numeric
    if(any( (is.numeric(value1) == FALSE) * (is.na(value1) == FALSE) )){
      stop(method, "wrong specification of \'", name1, "\' \n", 
           "\'", name1, "\' must be a numeric \n",        
           "is(", name1, ") : ", paste(is(value1), collapse = " "), "\n")
    }
    
    #### check duplicates
    if(refuse.duplicates == TRUE && any(duplicated(value1))){
      stop(method, "wrong specification of \'", name1, "\' \n", 
           "\'", name1, "\' contains duplicated values : \n",        
           paste(unique(value1[duplicated(value1)]), collapse = " "), "\n")
    }
    
    #### check min value
    if(!is.null(min) && any(stats::na.omit(value1) < min)){
      stop(method, "wrong specification of \'", name1, "\' \n", 
           "\'", name1, "\' must be bigger than ", min, " \n",        
           "invalid value(s) in ", name1, " : ", paste(value1[stats::na.omit(value1) < min], collapse = " "), "\n")
    }
    
    #### check max value
    if(!is.null(max) && any(stats::na.omit(value1) > max)){
      stop(method, "wrong specification of \'", name1, "\' \n", 
           "\'", name1, "\' must be smaller than ", max, " \n",        
           "invalid value(s) in ", name1, " : ", paste(value1[stats::na.omit(value1) > max], collapse = " "), "\n")
    }
    
    #### check valid values
    if(!is.null(validValues) && any(value1 %in% validValues == FALSE)){
      
      stop(method, "wrong specification of \'", name1, "\' \n", 
           "valid values for \'", name1, "\' : ", if(refuse.NULL == FALSE){"NULL"}, " \"", paste(validValues, collapse = "\" \""), "\" \n", 
           "refused value",if(sum(value1 %in% validValues == FALSE)>1){"s"}," for \'", name1, "\' : \"", paste(value1[value1 %in% validValues == FALSE], collapse = " "), "\"\n")
      
    }
  }
  
  return(invisible(TRUE))
}

#' @rdname validFCTs
#' @export
validPath <- function(value1, name1 = as.character(substitute(value1)), type,
                      method = NULL, addPP = TRUE, extension = NULL, checkFsep = FALSE){
  
  if(!is.null(method) && addPP){
    method <- paste0(method, ": ")
  }
  
  validCharacter(type, validLength = 1, validValues = c("file", "dir"))
  
  try_path <- switch(type,
                     file = file.exists(value1),
                     dir = dir.exists(value1)
  )
  
  if(try_path == FALSE){
    stop(method, "wrong specification of \'", name1, "\' \n", 
         "proposed ",type,": \"", value1, "\" is not valid \n", 
         "current path: ", getwd(), "\n")
  }
  
  
  if(type == "dir"){ 
    if(checkFsep == TRUE && substr(value1, start = nchar(value1), stop = nchar(value1)) != "/"){
    warning(method, "possible bad specification of \'", name1, "\' \n", 
            "\'", name1, "\' should end with a fsep (e.g. \"/\") \n", 
            "proposed ", type, " : ", value1, "\n")
    }
  }else if(type == "file" && !is.null(extension)){
    fileExtension <- tools::file_ext(value1) 
    if(fileExtension %in% extension == FALSE){
      stop(method, "\'", name1, "\' has not the expected extension \n", 
           "proposed extension: \"", fileExtension, "\" \n", 
           "expected extension: \"", paste(extension, collapse = "\" \""), "\"\n")
    }
  }
  
  return(invisible(TRUE))
}
