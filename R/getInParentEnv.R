### getInParentEnv.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt  5 2017 (10:48) 
## Version: 
## last-updated: okt  5 2017 (10:48) 
##           By: Brice Ozenne
##     Update #: 1
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#' @title Find object in the parent environments
#' 
#' @description Internal function
#' 
#' @param name character string containing the name of the object to get.
#' @param envir the environment from which to look for the object.
getInParentEnv <- function(name, envir){
  
  frames <- sys.status()
  all.frames <- sapply(1:length(frames$sys.frames), function(x){identical(parent.frame(x),globalenv())})
  index.parents <- which(all.frames==FALSE)
  n.parents <- length(index.parents)
  
  iParent <- 1
  res <- NULL
  while(iParent <= n.parents){ # iParent <- 1
    if(name %in% ls(envir = parent.frame(iParent))){
      res <- get(name, envir = parent.frame(iParent))
      iParent <- n.parents + 1
    }else{
      iParent <- iParent + 1     
    }
  }
  
  return(res)
}


#----------------------------------------------------------------------
### getInParentEnv.R ends here
