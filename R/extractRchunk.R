### extractRchunk.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 15 2017 (09:16) 
## Version: 
## Last-Updated: dec  2 2017 (12:27) 
##           By: Brice Ozenne
##     Update #: 93
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
#' @param file name of the file from which the R chunks should be extracted.
#' @param newfile name of the R file be created.
#' @param file.header string indicating headers.
#' @param file.start.SRC string indicating the start of a code block.
#' @param file.end.SRC string indicating the end of a code block.
#' @param rm.export.none should chunks of code with \code{:export none} be ignored?
#' @param index.chunk should the code block be numbered in the .R file?
#' @param overwrite Can the existing R file be overwritten?
#'
#' @export
extractRchunk <- function(file, newfile = NULL,
                          file.header = NULL, 
                          file.start.SRC = NULL,
                          file.end.SRC = NULL,
                          rm.export.none = TRUE,
                          index.chunk = TRUE, overwrite = FALSE){

### ** define tags
    type <- tolower(tools::file_ext(file))
    
    if(type == "org"){
        file.header <- "\\*"
        file.start.SRC <- "\\#\\+BEGIN_SRC R"
        file.end.SRC <- "\\#\\+END_SRC"

            if(length(grep(".org$",file))==0){
                stop("file must have a .org extension \n")
            }

            if(is.null(newfile)){
                newfile <- gsub(".org$",".R",file)
            }
            
        }else if(type == "rmd"){
            file.header <- "\\#"
            file.start.SRC <- "\\`\\`\\`\\{r"
            file.end.SRC <- "\\`\\`\\`$"

            if(length(grep(".rmd$",file))==0){
                stop("file must have a .rmd extension \n")
            }

            if(is.null(newfile)){
                newfile <- gsub(".rmd$",".R",file)
            }
            
        }else if(is.null(file.header)  || is.null(file.start.SRC) || is.null(file.end.SRC)){
            stop("if the file is not a org and rmd file the tags must be defined \n",
                 "please specify arguments \'file.header\', \'file.start.SRC\',  and \'file.end.SRC\' \n")
        }
        
    if(is.null(newfile)){
        stop("Please specify the name of the R file to be created using the argument \'newfile\' \n")
    }
    
### ** open file    
    if(file.exists(file)==FALSE){
        stop("file does not exist \n")
    }
    
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
        stop("Number of file.start.SRC does not match the number of file.end.SRC \n")
    }
    if(length(index.end.chunk) == 0){
        stop("No R code to extract \n")
        return(NULL)
    }

### ** group all
    df.extract <- rbind(data.frame(type = "header",
                                   index.start = index.header,
                                   index.stop = index.header,
                                   index = NA),
                        data.frame(type = "chunk",
                                   index.start = index.start.chunk,
                                   index.stop = index.end.chunk,
                                   index = 1:n.chunk))
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
### extractRchunk.R ends here
