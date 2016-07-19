._github_path <- "C:/Users/hpl802/Documents/GitHub/"


package.source <- function(name, path = ._github_path, Rcode = TRUE, Ccode = FALSE){
  
  validPath(path, type = "dir", method = "package.source")
  validPath(file.path(path, name), type = "dir", method = "package.source")
  
  if(Rcode){
    validPath(file.path(path, name, "R"), type = "dir", method = "package.source")
    R.utils::sourceDirectory(file.path(path, name, "R"), modifiedOnly = FALSE)
  }
  
  if(Ccode){
    validPath(file.path(path, name, "src"), type = "dir", method = "package.source")
    fileNames <- list.files(file.path(path, name, "src"))
    fileExts <- tools::file_path_sans_ext(file)
    indexC <- grep(".cpp|.c", x = tools::file_ext(lsFiles), fixed = FALSE)
    lapply(fileNames[indexC], Rcpp11::source.cpp)
  }
  
}
