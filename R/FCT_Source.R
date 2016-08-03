#' @title Find GitHub directory
#' 
#' @description Localise the github directory on linux or windows OS.
#'
#' @param user the name corresponding to the cession.
#' 
#' @export
dir.gitHub <- function(user = NULL){
  
  if(is.null(user)){
    if(Sys.info()["user"] != "unknown"){
      user <- Sys.info()["user"]
    }else if(Sys.info()["login"] != "unknown"){
      user <- Sys.info()["login"]
    }else stop("dir.gitHub: please specify the user \n")
  }
  
  if(Sys.info()["sysname"] == "Linux"){
    dir <- file.path("/home",user,"GitHub")
  }else if(Sys.info()["sysname"] == "Windows"){
    dir <- file.path("C:/Users",user,"Documents","GitHub")
  }else{
    stop("only implemented for linux and windows")
  }
  
  if(dir.exists(dir)){
    return(dir) 
  }else{
    stop("dir.gitHub: no GitHub directory found")
  }
  
}

#' @title Source a package directory
#' 
#' @description Source all the R and Cpp file contain in a package
#' 
#' @param name the name of the package
#' @param path the path to the directory containing the package
#' @param Rcode should all the .R be sourced
#' @param Ccode should all the .cpp file be source (using Rcpp::sourceCpp)
#'
#' @seealso \code{\link{dir.gitHub}}
#'
#' @export
package.source <- function(name, path = dir.gitHub(), Rcode = TRUE, Ccode = FALSE){
  
  validPath(path, type = "dir", method = "package.source")
  validPath(file.path(path, name), type = "dir", method = "package.source")
  
  if(Rcode){
    validPath(file.path(path, name, "R"), type = "dir", method = "package.source")
    R.utils::sourceDirectory(file.path(path, name, "R"), modifiedOnly = FALSE, envir = globalenv())
  }
  
  if(Ccode){
    validPath(file.path(path, name, "src"), type = "dir", method = "package.source")
	fileNames <- list.files(file.path(path, name, "src"))
	fileExts <- tools::file_path_sans_ext(fileNames)
	indexC <- grep("cpp", x = tools::file_ext(fileNames), 
					fixed = FALSE)
	lapply(file.path(path,name,"src",setdiff(fileNames[indexC],"RcppExports.cpp")), 
		   Rcpp::sourceCpp, 
		   rebuild = TRUE)

  }
  
}
