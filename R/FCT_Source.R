#' @title Source a package directory
#' 
#' @description Source all the R and Cpp file contain in a package
#' 
#' @param name the name of the package
#' @param path the path to the directory containing the package
#' @param Rpackage should the related R package be loaded
#' @param Rcode should all the .R be sourced
#' @param RorderDescription should the R files be sourced in the order indicate by collate
#' @param attach source the .onAttach function if it is present in the current environment
#' @param Ccode should all the .cpp file be source (using Rcpp::sourceCpp)
#' @param rebuild Force a rebuild of the shared library (from Rcpp:::sourceCpp).
#' @param warning should a warning be displayed if some of the R files are not sourced
#' 
#' @seealso \code{\link{path_gitHub}}
#'
#' @examples 
#' package.source("butils")
#'
#' @export
package.source <- function(name, path = path_gitHub(), 
                           Rpackage = TRUE,
                           Rcode = TRUE, RorderDescription = TRUE, attach = TRUE,
                           Ccode = FALSE, rebuild = FALSE,
                           warning = TRUE){
  
  validPath(path, type = "dir", method = "package.source")
  validPath(file.path(path, name), type = "dir", method = "package.source")
  if(Rpackage){
    
    if(identical(Rpackage,TRUE)){Rpackage <- c("Imports","Depends","Suggests")}
    
    validCharacter(Rpackage, validLength = NULL, validValues = c("Imports","Depends","Suggests"))
    
    packageToLoad <- NULL
    if("Imports" %in% Rpackage){
      packageToLoad <- c(packageToLoad,
                         read_description(name, path = path, field = "Imports", rmComma = TRUE, rmBlanck = TRUE)
      )
    }
    if("Depends" %in% Rpackage){
      packageToLoad <- c(packageToLoad,
                         read_description(name, path = path, field = "Depends", rmComma = TRUE, rmBlanck = TRUE)
      )
    }
    if("Suggests" %in% Rpackage){
      packageToLoad <- c(packageToLoad,
                         read_description(name, path = path, field = "Suggests", rmComma = TRUE, rmBlanck = TRUE)
      )
    }
    
    ## remove version in parenthesis
    packageToLoad <- gsub( " *\\(.*?\\) *", "", packageToLoad)
    
    packageToLoad <- setdiff(packageToLoad, "R")
    
    lapply(packageToLoad, require, character.only = TRUE)
  }
  
  if(Rcode){
    validPath(file.path(path, name, "R"), type = "dir", method = "package.source")
    
    ## find files
    fileNames <- setdiff(list.files(file.path(path, name, "R")),
                         "RcppExports.R")
    fileExts <- tools::file_path_sans_ext(fileNames)
    indexC <- grep("R", x = tools::file_ext(fileNames), 
                   fixed = FALSE)
    fileNames <- fileNames[indexC]
      
    if(RorderDescription){ ## reorder according DESCRIPTION
      
      if(file.exists(file.path(path,name,"DESCRIPTION")) == FALSE){
        warning("package.source: no DESCRIPTION file founded \n",
                "set \'RorderDescription\' to FALSE to source all the files that are present in the directory R \n")  
      }
      
      file.description <- readLines(file.path(path,name,"DESCRIPTION"))
      indexLine <- grep("Collate:",file.description)+1
      indexLine_end <-  min(grep(":",file.description,fixed=TRUE)[grep(":",file.description,fixed=TRUE)>indexLine])-1
      filesR.description <- file.description[indexLine:indexLine_end]
      filesR.description <- gsub("[[:blank:]]|'", "", filesR.description)
      
      test.missing <- is.na(match(fileNames, filesR.description))
      if(warning && any(test.missing)){
        warning("package.source: did not find files: ",paste(fileNames[which(test.missing)], collapse = " ")," in DESCRIPTION \n",
                "set \'RorderDescription\' to FALSE to source all the files that are present in the directory R \n")  
      }
      fileNames <- filesR.description[filesR.description %in% fileNames]
    }
    
    .onAttach_save <- NULL
    
    ## SOURCE
    lapply(file.path(path,name,"R",fileNames), source)
    
    ## mimic .onload
    if(attach && exists(".onAttach") && !is.null(.onAttach)){
      .onAttach()
    }
  }
  
  if(Ccode){
    validPath(file.path(path, name, "src"), type = "dir", method = "package.source")
    fileNames <- list.files(file.path(path, name, "src"))
    fileExts <- tools::file_path_sans_ext(fileNames)
    indexC <- grep("cpp", x = tools::file_ext(fileNames), 
                   fixed = FALSE)
    lapply(file.path(path,name,"src",setdiff(fileNames[indexC],"RcppExports.cpp")), 
           Rcpp::sourceCpp, 
           rebuild = rebuild)
    
  }
  
  return(invisible(TRUE))
  
}



