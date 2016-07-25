#' @title localise the GitHub directory
#'
#' @export
dir.gitHub <- function(user = NULL){
  
  if(Sys.info()["sysname"] == "linux"){
    
  }else if(Sys.info()["sysname"] == "Windows"){
    
    if(is.null(user)){user <- Sys.info()["login"]}
    dir <- file.path("C:/Users",user,"Documents","GitHub")
    
    if(dir.exists(dir)){
      return(dir) 
    }else{
      stop("dir.gitHub: no GitHub directory found")
    }
    
  }else{
    
    stop("only implemented for linux and windows")
    
  }
  
}

#' @title Source a package directory
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
