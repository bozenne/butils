profileCode <- function(FUN, filename = NULL, type = "Rprof", plot = TRUE){
  if(is.null(filename)){
    filename <- tempfile()
  }
  
  if(type == "Rprof"){
    
    Rprof(filename, line.profiling = TRUE)
    FUN()
    Rprof(append = FALSE)
    print(sProfile <- summaryRprof(filename))
    resProfile <- proftools::readProfileData(filename)
    unlink(filename)
    
    if(plot){
      proftools::plotProfileCallGraph(resProfile, score = "total", nodeSizeScore = "none")
    }
    
  }else if(type == "lineprof"){
    sProfile <- lineprof::lineprof(fctTest())
    if(plot){lineprof::shine(sProfile)}
  }
  
  
  
  
  return(invisible(sProfile))
}
