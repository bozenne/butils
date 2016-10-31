#' @title Load description file
#' @description Extract the description file of a package
#' 
#' @param package the name of the package
#' @param path the path to the directory containing the package
#' @param field should a specific field of the description file be extracted
#' @param rmBlanck should the initial blanck be removed?
#' @param rmComma should any comma be removed?
#'
#' @examples 
#' \dontrun{
#' read_description("butils", path = path_gitHub())
#' read_description("butils", path = path_gitHub(), field = "Collate")
#' read_description("butils", path = path_gitHub(), field = "Collate", rmBlanck = FALSE)
#' 
#' read_description("butils", path = path_gitHub(), field = "Version")
#' read_description("butils", path = path_gitHub(), field = "Imports")
#' }
read_description <- function(package, path, field = NULL,
                             rmBlanck = TRUE, rmComma = TRUE){
  
  validPath(path, type = "dir", method = "read_description")
  validPath(file.path(path, package), type = "dir", method = "read_description")
  
  file.description <- readLines(file.path(path,package,"DESCRIPTION"))
  
  if(is.null(field)){return(file.description)}
  
  # select field
  if(length(grep(field,file.description))==0){return(NULL)}
  indexLine <- grep(field,file.description)
  indexLine_end <-  min(grep(":",file.description,fixed=TRUE)[grep(":",file.description,fixed=TRUE)>indexLine])-1
  file.description <- file.description[seq(indexLine, indexLine_end)]
  
  # remove initial blanck
  if(rmBlanck){
    file.description[1] <- gsub(paste0(field,":"),"",file.description[1]) # remove field name
     
    file.description <- trimws(file.description, which = "left") # remove blancks
    file.description <- file.description[sapply(file.description, nchar)>0] # remove first line
  }
  # remove comma
  if(rmComma){
    file.description <- gsub(",", "", x = file.description)
  }
  
  # export
  return(file.description)
  
}

#' @title Update description file
#' @description Overwrite the description file of a package
#' 
#' @inheritParams read_description
#' @param newfile the content of the new description file
#' @param trace should the execution of the function be traced
#'
#' @examples 
#' \dontrun{
#' file <- read_description("butils", path = path_gitHub())
#' write_description("butils", path = path_gitHub(), file)
#' }
write_description <- function(package, path, newfile, trace = TRUE){
  
  validPath(path, type = "dir", method = "write_description")
  validPath(file.path(path, package), type = "dir", method = "write_description")
  validPath(file.path(path, package, "DESCRIPTION"), type = "file", method = "write_description")
  
  ## write
  con <- file(file.path(path,package,"DESCRIPTION")) 
  writeLines(text = newfile, con = con) 
  close(con)
  
  return(invisible(TRUE))
}


#' @title write cpllate field
#' @description write thecpllate field in the DESCRIPTION file of a package
#' 
#' @inheritParams read_description
#' @param trace should the execution of the function be traced
#' 
#' @examples 
#' \dontrun{
#' write_collate("butils", path = path_gitHub())
#' }
write_collate <- function(package, path, trace = TRUE){
  
  ## read description file
  file.description <- read_description(package, path = path)
  
  ## count white space
  line1 <- read_description(package, path = path, field = "Collate", rmBlanck = FALSE)[2] # 2 because the first line is the field name
  line0 <- read_description(package, path = path, field = "Collate", rmBlanck = TRUE)[1]
  space <- stringr::str_pad("", width=nchar(line1)-nchar(line0), side="right")
  
  ## list all R files
  Rfiles <- list.files(file.path(path,package,"R"), pattern = c(".R$",".r$"))
  
  ## new collapse
  TOadd <- c("Collate:",paste0(space,"'",Rfiles,"'"))
  
  indexLine <- grep("Collate:",file.description)
  test.change  <- FALSE 
  
  if(length(indexLine) == 0){
    if(trace){cat(">> add \'collate\' field to the DESCRIPTION file")}
    file.description <- c(file.description, TOadd)
    test.change <- TRUE
  }else{
    testPP <- grep(":",file.description,fixed=TRUE)>indexLine
    if(any(testPP)){
      indexLine_end <- min(grep(":",file.description,fixed=TRUE)[testPP])-1
    }else{
      indexLine_end <-  length(file.description)
    }
    
    if(!identical(TOadd,file.description[indexLine:indexLine_end])){
      if(trace){cat(">> update \'collate\' field in the DESCRIPTION file")}
      file.description <- c(file.description[1:(indexLine-1)],TOadd,file.description[(indexLine_end+1):length(file.description)])
      test.change <- TRUE
    }
  }
  
  if(test.change){
    write_description(package, path = path, newfile = file.description)
    if(trace){cat("\n")}
  }
  return(invisible(TRUE))
}

#' @title write date field
#' @description write the date field in the DESCRIPTION file of a package
#' 
#' @inheritParams read_description
#' @param trace should the execution of the function be traced
#' 
#' @examples 
#' \dontrun{
#' write_date("butils", path = path_gitHub())
#' }
write_date <- function(package, path, trace = TRUE){
  
  ## read description file
  file.description <- read_description(package, path = path)
  indexLine <- grep("Date:",file.description)
  
  newDate <- paste0("Date: ",format(Sys.time(), "%Y-%M-%d"))
  
  if(length(indexLine) == 0){
    warning("writeDate_package: \'Date\' field is missing in DESCRIPTION \n")
  }else if(file.description[indexLine] != newDate){
    if(trace){cat(">> update \'Date\' field in the DESCRIPTION file")}
    file.description[indexLine] <- paste0("Date: ",format(Sys.time(), "%Y-%M-%d"))
    write_description(package, path = path, newfile = file.description)
    if(trace){cat("\n")}
  }
  
  return(invisible(TRUE))
}