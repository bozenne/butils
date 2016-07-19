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
#' @rdname validFCTs
#' @export
validCharacter <- function(value, name = as.character(substitute(value)), validLength, 
                           validValues = "character", refuse.NULL = TRUE, refuse.duplicates = FALSE, 
                           method = NULL, addPP = TRUE){
  
  if(!is.null(method) && addPP){
    method <- paste0(method, ": ")
  }
  
  if(is.null(value)){
    
    if(refuse.NULL == TRUE){
      stop(method, "wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must not be NULL \n")
    }
    
  }else{
    
    #### check size
    n.value <- length(value)
    
    if(!is.null(validLength) && n.value %in% validLength == FALSE){
      stop(method, "wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must have length ", paste(validLength, collapse = " or "), "  \n", 
           "length(", name, ") : ", n.value, "\n")
    }
    
    #### check duplicates
    if(refuse.duplicates == TRUE && any(duplicated(value))){
      stop(method, "wrong specification of \'", name, "\' \n", 
           "\'", name, "\' contains duplicated values : \n",        
           "\"",paste(unique(value[duplicated(value)]), collapse = "\" \""), "\" \n")
    }
    
    #### check values
    if(identical(validValues,"character")){
      
      if(any(is.character(value) == FALSE)){
        stop(method, "wrong specification of \'", name, "\' \n", 
             "\'", name, "\' must be a ", if(n.value == 1){"character"}else{"vector of characters"}," \n", 
             "is(", name, ") : ", paste(is(value), collapse = " "), "\n")
      }
      
    } else if(identical(validValues,"character_or_logical")){
      
      if(any( (is.character(value) == FALSE) * (is.logical(value) == FALSE) > 0 )){
        stop(method, "wrong specification of \'", name, "\' \n", 
             "\'", name, "\' must be a ", if(n.value == 1){"character or logical"}else{"vector of characters or logicals"}," \n", 
             "is(", name, ") : ", paste(is(value), collapse = " "), "\n")
      }
      
    } else if(!is.null(validValues) && any(value %in% validValues == FALSE)){
      
      stop(method, "wrong specification of \'", name, "\' \n", 
           "valid values for \'", name, "\' : ", if(refuse.NULL == FALSE){"NULL"}, " \"", paste(validValues, collapse = "\" \""), "\" \n", 
           "refused value",if(sum(value %in% validValues == FALSE)>1){"s"}," for \'", name, "\' : \"", paste(value[value %in% validValues == FALSE], collapse = "\" \""), "\"\n")
      
    }
    
  }
  
  return(invisible(TRUE))
  
}

#' @rdname validFCTs
#' @export
validClass <- function(value, name = as.character(substitute(value)), validClass, 
                       superClasses = TRUE, method = NULL, addPP = TRUE){
  
  if(!is.null(method) && addPP){
    method <- paste0(method, ": ")
  }
  
  if(superClasses == TRUE){
    
    if( all(is(value) %in% validClass == FALSE) ){
      stop(method, "wrong specification of \'", name, "\' \n", 
           "superclass of \'", name, "\' must be one of the following \"", paste(validClass,collapse="\" \""), "\"  \n", 
           "proposed superclass : \"", paste(is(value),collapse="\" \""), "\" \n")
    }  
    
  }else{
 
    if( class(value) %in% validClass == FALSE){
      stop(method, "wrong specification of \'", name, "\' \n", 
           "class of \'", name, "\' must be \"", paste(validClass,collapse="\" \""),"\"  \n", 
           "proposed class : ", class(value)[[1]], "\n")
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
validInteger <- function(value, name = as.character(substitute(value)), validLength, 
                         validValues = NULL, min = NULL, max = NULL, 
                         refuse.NA = TRUE, refuse.NULL = TRUE, refuse.duplicates = FALSE, 
                         method = NULL, addPP = TRUE){
  
  if(!is.null(method) && addPP){
    method <- paste0(method, ": ")
  }
  
  validNumeric(value = value, name = name, validLength = validLength, min = min, max = max, 
               refuse.NA = refuse.NA, refuse.NULL = refuse.NULL, refuse.duplicates = refuse.duplicates, method = method)
  
  #### check integer
  if(!is.null(value) && any(value %% 1 > 0)){
    stop(method, "wrong specification of \'", name, "\' \n", 
         "\'", name, "\' must contain integers not doubles \n",        
         "invalid value(s) in ", name, " : ", paste(value[value %% 1 > 0], collapse = " "), "\n")
  }
  
  return(invisible(TRUE))
}

#' @rdname validFCTs
#' @export
validLogical <- function(value, name = as.character(substitute(value)), validLength, 
                         refuse.NULL = TRUE, refuse.NA = TRUE, 
                         method = NULL, addPP = TRUE){
  
  if(!is.null(method) && addPP){
    method <- paste0(method, ": ")
  }
  
  if(is.null(value)){
    
    #### NULL
    if(refuse.NULL == TRUE){
      stop(method, "wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must be logical ",if(refuse.NA == FALSE){"or NA"}," and not NULL \n")
    }
    
  }else{ 
    
    #### Size
    if(!is.null(validLength) && length(value) %in% validLength == FALSE){
      stop(method, "wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must have length ", paste(validLength, collapse = " or "), "  \n", 
           "length(", name, ") : ", length(value), "\n")
    } 
    
    #### Type
    if(any(is.logical(value) == FALSE)){
      stop(method, "wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must be ", if(refuse.NULL == FALSE){"NULL or "}, if(refuse.NA == FALSE){"NA or "},"TRUE or FALSE \n",        
           "is(", name, ") : ", paste(is(value), collapse = " "), "\n")
    }
    
    if(refuse.NA == TRUE && any(is.na(value)) ){
      stop(method, "wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must be logical ",if(refuse.NULL == FALSE){"or NULL"}," and not NA \n")
    }
    
  }
  
  return(invisible(TRUE))
}

#' @rdname validFCTs
#' @export
validNames <- function(value, name = as.character(substitute(value)), refuse.NULL = TRUE,
                       validLength = NULL, validValues = NULL, requiredValues = NULL, forbiddenValues = NULL,
                       method = NULL, addPP = TRUE){
  
  if(!is.null(method) && addPP){
    method <- paste0(method, ": ")
  }
  
  ## type
  if(is.matrix(value)){
    value <- colnames(value)
  }
  
  if(is.data.table(value) || is.data.frame(value) || is.list(value)){
    value <- names(value)
  }
  
  ## tests
  if(is.null(value)){
    
    if(refuse.NULL == TRUE){
    stop(method, "wrong specification of \'", name, "\' \n", 
         "names of \'", name, "\' must not be NULL \n")
    }
    
  }else{
    
    #### check size
    n.value <- length(value)
    
    if(!is.null(validLength) && n.value %in% validLength == FALSE){
      stop(method, "wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must have ", paste(validLength, collapse = " or ")," names  \n", 
           "length(names(", name, ")) : ", n.value, "\n")
    }
    
    #### check content
    
    if(!is.null(requiredValues) && any(requiredValues %in% value == FALSE)){
      
      stop(method, "wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must contains specific names \n",
           "missing names : \"",paste(requiredValues[requiredValues %in% value == FALSE], collapse = "\" \""),"\" \n", 
           "proposed names : \"", paste(value, collapse = "\" \""), "\"\n")  
      
    }
    
    if(!is.null(validValues) && any(value %in% validValues == FALSE)){
      
      stop(method, "wrong specification of \'", name, "\' \n", 
           "valid names for \'", name, "\' : \"",paste(validValues, collapse = "\" \""),"\" \n", 
           "refused names : \"", paste(value[value %in% validValues == FALSE], collapse = " "), "\"\n")  
      
    }
    
    if(!is.null(forbiddenValues) && any(value %in% forbiddenValues)){
      
      stop(method, "wrong specification of \'", name, "\' \n", 
           "forbidden names for \'", name, "\' : \"",paste(forbiddenValues, collapse = "\" \""),"\" \n", 
           "refused names : \"", paste(value[value %in% forbiddenValues], collapse = " "), "\"\n")  
      
    }
    
    if(any(duplicated(value))){
      stop(method, "wrong specification of \'", name, "\' \n", 
           name, " must not contain duplicated names \n", 
           "duplicated names : \"", paste(value[duplicated(value)], collapse = " "), "\"\n")  
    }
    
  }
  
  return(invisible(TRUE))
  
}
#' @rdname validFCTs
#' @export
validNumeric <- function(value, name = as.character(substitute(value)), validLength,
                         validValues = NULL , min = NULL, max = NULL,
                         refuse.NA = TRUE, refuse.NULL = TRUE, refuse.duplicates = FALSE, 
                         method = NULL, addPP = TRUE){
  
  if(!is.null(method) && addPP){
    method <- paste0(method, ": ")
  }
  
  if(is.null(value)){
    
    if(refuse.NULL == TRUE){
      stop(method, "wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must not be NULL \n")
    }
    
  }else{
    
    #### check length
    if(!is.null(validLength) && length(value) %in% validLength == FALSE){
      stop(method, "wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must have length ", paste(validLength, collapse = " or "), "  \n", 
           "length(", name, ") : ", length(value), "\n")
    }
    
    #### check NA
    if(refuse.NA == TRUE && any(is.na(value))){
      stop(method, "wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must not contain NA \n", 
           "index of NA values : ", which(paste(is.na(value), collapse = " ")), "\n")
    }
    
    #### check numeric
    if(any( (is.numeric(value) == FALSE) * (is.na(value) == FALSE) )){
      stop(method, "wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must be a numeric \n",        
           "is(", name, ") : ", paste(is(value), collapse = " "), "\n")
    }
    
    #### check duplicates
    if(refuse.duplicates == TRUE && any(duplicated(value))){
      stop(method, "wrong specification of \'", name, "\' \n", 
           "\'", name, "\' contains duplicated values : \n",        
           paste(unique(value[duplicated(value)]), collapse = " "), "\n")
    }
    
    #### check min value
    if(!is.null(min) && any(stats::na.omit(value) < min)){
      stop(method, "wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must be bigger than ", min, " \n",        
           "invalid value(s) in ", name, " : ", paste(value[stats::na.omit(value) < min], collapse = " "), "\n")
    }
    
    #### check max value
    if(!is.null(max) && any(stats::na.omit(value) > max)){
      stop(method, "wrong specification of \'", name, "\' \n", 
           "\'", name, "\' must be smaller than ", max, " \n",        
           "invalid value(s) in ", name, " : ", paste(value[stats::na.omit(value) > max], collapse = " "), "\n")
    }
    
    #### check valid values
    if(!is.null(validValues) && any(value %in% validValues == FALSE)){
      
      stop(method, "wrong specification of \'", name, "\' \n", 
           "valid values for \'", name, "\' : ", if(refuse.NULL == FALSE){"NULL"}, " \"", paste(validValues, collapse = "\" \""), "\" \n", 
           "refused value",if(sum(value %in% validValues == FALSE)>1){"s"}," for \'", name, "\' : \"", paste(value[value %in% validValues == FALSE], collapse = " "), "\"\n")
      
    }
  }
  
  return(invisible(TRUE))
}

#' @rdname validFCTs
#' @export
validPath <- function(value, name = as.character(substitute(value)), type,
                      method = NULL, addPP = TRUE, extension = NULL, checkFsep = FALSE){
  
  if(!is.null(method) && addPP){
    method <- paste0(method, ": ")
  }
  
  validCharacter(type, validLength = 1, validValues = c("file", "dir"))
  
  try_path <- switch(type,
                     file = file.exists(value),
                     dir = dir.exists(value)
  )
  
  if(try_path == FALSE){
    stop(method, "wrong specification of \'", name, "\' \n", 
         "proposed ",type,": \"", value, "\" is not valid \n", 
         "current path: ", getwd(), "\n")
  }
  
  
  if(type == "dir"){ 
    if(checkFsep == TRUE && substr(value, start = nchar(value), stop = nchar(value)) != "/"){
    warning(method, "possible bad specification of \'", name, "\' \n", 
            "\'", name, "\' should end with a fsep (e.g. \"/\") \n", 
            "proposed ", type, " : ", value, "\n")
    }
  }else if(type == "file" && !is.null(extension)){
    fileExtension <- tools::file_ext(value) 
    if(fileExtension %in% extension == FALSE){
      stop(method, "\'", name, "\' has not the expected extension \n", 
           "proposed extension: \"", fileExtension, "\" \n", 
           "expected extension: \"", paste(extension, collapse = "\" \""), "\"\n")
    }
  }
  
  return(invisible(TRUE))
}
