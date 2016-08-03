#' @title Profile a function
#' 
#' @description This function is a convenient to profile the execution of a function
#' 
#' @param FUN the function to be profiled. Must have no argument
#' @param filename the name of the file where to save the profiling. I
#' @param type the function used for profiling: can be \code{"Rprof"} or \code{"lineprof"}.
#' @param plot should the result of the profiling be graphically displayed. Only relevant if \code{type="Rprof"}.
#' 
#' @details 
#' If \code{filename} is \code{NULL} a temporary file is created and is removed at the end of the execution of the function.
#'   
#' @return invisible output of the profiling method.
#' 
#' @keywords function profile
#' @export
profileCode <- function(FUN, filename = NULL, type = "Rprof", plot = TRUE){
  if(is.null(filename)){
    filename <- tempfile()
    on.exit(unlink(filename))
  }
  
  if(type == "Rprof"){
    Rprof(filename, line.profiling = TRUE)
    FUN()
    Rprof(append = FALSE)
    print(sProfile <- summaryRprof(filename))
    resProfile <- proftools::readProfileData(filename)
    
    if(plot){
      proftools::plotProfileCallGraph(resProfile, score = "total", nodeSizeScore = "none")
    }
    
  }else if(type == "lineprof"){
    sProfile <- lineprof::lineprof(fctTest())
    if(plot){lineprof::shine(sProfile)}
  }
  
  
  
  
  return(invisible(sProfile))
}
