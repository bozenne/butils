### extractSRCorgmode.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 15 2017 (09:16) 
## Version: 
## Last-Updated: nov 24 2017 (09:18) 
##           By: Brice Ozenne
##     Update #: 59
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

#' @title Extract code chunks from an .org file
#' @description Collect code chunks from an .org file and put them into an .R file.
#'
#' @param file name of the org file
#' @param file.header string indicating headers
#' @param file.start.SRC string indicating the start of a code block
#' @param file.end.SRC string indicating the end of a code block
#' @param rm.export.none should chunks of code with \code{:export none} be ignored?
#' @param index.chunk should the code block be numbered in the .R file?
#' @param overwrite Can the existing R file be overwritten?
#' 
#' @export
extractSRCorgmode <- function(file,
                              file.header = "\\*", 
                              file.start.SRC = "\\#\\+BEGIN_SRC R",
                              file.end.SRC = "\\#\\+END_SRC",
                              rm.export.none = TRUE,
                              index.chunk = TRUE, overwrite = FALSE){


### ** open file
    if(length(grep(".org$",file))==0){
        stop("file must have a .org extension \n")
    }
    if(file.exists(file)==FALSE){
        stop("file does not exist \n")
    }
    newfile <- gsub(".org$",".R",file)
    if(file.exists(newfile) && (overwrite == FALSE)){
        stop("corresponding R file already exists \n",
             "set argument \'overwrite\' to TRUE to overwrite it \n")
    }
    
    con <- file(file, "rb") 
    file.line <- readLines(con) 
    close(con)

### ** find headers (if any)
    index.header <- grep(paste0("^",file.header),file.line, value = FALSE)
    n.header <- length(index.header)
        
### ** find chunk (if any)
    index.start.chunk <- grep(paste0("^",file.start.SRC),
                              x = file.line, value = FALSE)
    index.end.chunk <- grep(paste0("^",file.end.SRC),
                            x = file.line, value = FALSE)
    
    n.chunk <- length(index.start.chunk)
    if(length(index.end.chunk) != n.chunk){
        stop("Number of file.start.SRC does not match the number of file.end.src \n")
    }

### ** group all
    df.extract <- rbind(data.frame(type = "header", index.start = index.header, index.stop = index.header, index = NA),
                        data.frame(type = "chunk", index.start = index.start.chunk, index.stop = index.end.chunk, index = 1:n.chunk))
    df.extract <- df.extract[order(df.extract$index.start),]
    n.extract <- NROW(df.extract)
    
### ** create string
    file.content <- NULL
    
    for(iE in 1:n.extract){ # iE <- 1
        iType <- as.character(df.extract[iE,"type"])
        iStart <- df.extract[iE,"index.start"]
        iEnd <- df.extract[iE,"index.stop"]
        if(iType == "chunk"){
            test.export.none <- length(grep(":exports none",file.line[iStart],fixed = TRUE))
            if(rm.export.none && test.export.none==1){
                    next
            }
            iStart <- iStart + 1
            iEnd <- iEnd - 1
        }
        
        add.before <- switch(iType,
                             "header"=NULL,
                             "chunk"=if(index.chunk){paste0("## chunk ",df.extract[iE,"index"],"")}else{NULL})

        prefix.add <- switch(iType,
                             "header"="## ",
                             "chunk"=NULL)
        toAdd <- paste0(prefix.add,file.line[iStart:iEnd])
        
        add.after <- switch(iType,
                            "header"=if(iE<n.extract&&df.extract[iE+1,"type"]=="chunk"){""}else{NULL},
                            "chunk"="")
        
        file.content <- c(file.content,add.before,toAdd,add.after)
    }

### ** write file
    con <- file(newfile) 
    writeLines(text = file.content, con = con) 
    close(con)

### ** export
    return(invisible(file.content))
    
    
}



##----------------------------------------------------------------------
### extractSRCorgmode.R ends here
