#' @title Launching batch 
#' 
#' @description Convenient way to launch an R file in batch mode.
#' 
#' @param file the name of the file. \emph{character}.
#' @param path the path to the directory containing the file. \emph{character}.
#' @param dirBatch If no \code{NULL}, the name of a directory where to store the .lis and .Rout. \emph{character}.  
#' @param operator the command used to launch the batch. \emph{character}.
#' @param add.Rout Should a file .Rout be generated. \emph{character}.  
#' @param add.lis Should a file .liss be generated.\emph{character}.  
#' @param add.options Additional options to be used to launch the batch. \emph{character}.  
#' @param rm.newfile Should the copy of the original file be removed after completion of the function? \emph{character}.  
#' 
#' @details 
#' If the directory specified for argument \code{dirBatch} does not exist, it will be automatically generated (provided it corresponds to a valid path).
#' Before launching the batch, first a copy of the file is made. It is named with the prefix \code{launchBatch-} and, if any, will automatically overwrite the existing file.
#' Then the command for generating the .lis is added to the file and the batch is launched.
#'   
#' @return the output of the shell when lauching the batch. \emph{character}.  
#' 
#' @keywords function batch
#' @export
launchBatch <- function(file, path = ".", dirBatch = NULL,  
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
