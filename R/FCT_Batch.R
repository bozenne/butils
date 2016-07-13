#' @title Launching batch 
#'
launchBatch <- function(path = ".", file, dirBatch = NULL, txt.options, 
                        operator = "start R CMD BATCH", add.Rout = TRUE, add.lis = TRUE, add.options = "",
                        rm.newfile = FALSE){
  
  validPath(path, type = "dir", method = "launchBatch")
  validPath(file.path(path, file), type = "file", method = "launchBatch", extension = c("r","R"))
  
  #### analyse file path
  fileName <- tools:::file_path_sans_ext(file)
  fileExtension <- tools:::file_ext(file)
  
  #### directory for the batch file
  if(!is.null(dirBatch)){
    pathDir <- file.path(path,dirBatch)
    if(!dir.exists(pathDir)){
      dir.create(pathDir)
    }
  }else{
    pathDir <- path
  }
  
  #### temporary copy
  newfile <- file.path(pathDir,paste0("launchBatch-",fileName,".",fileExtension))
  
  file.copy(from = file.path(path,file), to = newfile, overwrite = TRUE)
  if(rm.newfile){on.exit(file.remove(newfile))}
  
  #### add liss
  if(add.lis){
    txt <- readLines(newfile, warn = FALSE)
    fileLis <- file.path(pathDir,paste0("launchBatch-",fileName,".lis"))
     
    newtxt <- c(paste0("sink(\"",fileLis,"\")"), 
                txt,
                "sink()")
    writeLines(newtxt, con = newfile)
  }
  
  #### command
  command <- paste0(operator, " \"", newfile,"\"")
  if(add.Rout){
    command <- paste(command, paste0("\"",newfile,"out\""))
  }
  command <- paste0(command, add.options)
  
  
  
  output <- shell(command)
  return(invisible(output))
}
