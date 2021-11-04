### extractRchunk.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 15 2017 (09:16) 
## Version: 
## Last-Updated: Nov  4 2021 (22:39) 
##           By: Brice Ozenne
##     Update #: 118
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
#' @param rm.noexport should header with \code{:noexport:} be ignored?
#' @param index.chunk should the code block be numbered in the .R file?
#' @param overwrite Can the existing R file be overwritten?
#'
#' @export
extractRchunk <- function(file, newfile = NULL,
                          file.header = NULL, 
                          file.start.SRC = NULL,
                          file.end.SRC = NULL,
                          rm.export.none = TRUE,
                          rm.noexport = TRUE,
                          index.chunk = TRUE,
                          overwrite = FALSE){

### ** define tags
    type <- tolower(tools::file_ext(file))
    
    if(type == "org"){
        file.header <- "\\*"
        file.start.SRC <- "\\#\\+BEGIN_SRC R"
        file.end.SRC <- "\\#\\+END_SRC"

        if(identical(rm.noexport,TRUE)){
            rm.noexport <- ":noexport:"
        }else if(identical(rm.noexport,FALSE)){
            rm.noexport <- NULL
        }
        
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

        if(identical(rm.noexport,TRUE)){
            rm.noexport <- NULL
        }else if(identical(rm.noexport,FALSE)){
            rm.noexport <- NULL
        }
        
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
        if(length(index.start.chunk)>length(index.end.chunk)){
            ## when the next chunk start before the end of the previous one
            index.pb <- which(index.start.chunk[-1]<c(index.end.chunk,rep(0,length(index.start.chunk)-length(index.end.chunk)-1)))
            if(length(index.pb)>0){
                stop("Number of file.start.SRC exceed the number of file.end.SRC. \n",
                     "Possible issue between line: ",index.start.chunk[-1][index.pb[1]]," and ",index.end.chunk[index.pb[1]],"\n")
            }else{
                stop("Number of file.start.SRC exceed the number of file.end.SRC. \n")
            }
        }else{
            stop("Number of file.end.SRC exceed the number of file.start.SRC. \n")
        }
        
    }
    if(length(index.end.chunk) == 0){
        stop("No R code to extract \n")
        return(NULL)
    }

### ** group all
    dt.extract <- rbind(data.table(type = "header",
                                   index.start = index.header,
                                   index.stop = index.header,
                                   index = NA),
                        data.table(type = "chunk",
                                   index.start = index.start.chunk,
                                   index.stop = index.end.chunk,
                                   index = 1:n.chunk))
    setkeyv(dt.extract, "index.start")
    n.extract <- NROW(dt.extract)

### ** level of the headings
    indexDT.header <- which(dt.extract$type == "header")
    index.header <- dt.extract[indexDT.header]$index.start
    vec.star <- unlist(lapply(strsplit(file.line[index.header], split = "* "),"[[",1))
    dt.extract[, c("level") := as.integer(NA)]
    dt.extract[indexDT.header, c("level") := nchar(vec.star)]

### ** end of the headings
    dt.extract[,c("index.stopHeader") := as.integer(NA)]
    dt.extract[indexDT.header, c("index.stopHeader") := c(.SD$index.start[-1], length(file.line)), by = "level"]

### ** remove non exported sections
    dt.extract[, c("export") := TRUE]
    if(!is.null(rm.noexport)){
        indexNoexport.header <- grep(":noexport:",file.line[index.header])
        if(length(indexNoexport.header)>0){
            
            dt.extract[indexDT.header[indexNoexport.header], c("export") := FALSE ]
            
            ls.index <- mapply(dt.extract[dt.extract$export==FALSE, .SD$index.start],
                               dt.extract[dt.extract$export==FALSE, .SD$index.stopHeader],
                               FUN=seq, SIMPLIFY = FALSE)
            vec.noexport <- sort(unique(unlist(ls.index)))

            dt.extract[dt.extract$index.start %in% vec.noexport, c("export") := FALSE]
        }
       
    }
    
### ** create string
    file.content <- NULL
    
    for(iE in 1:n.extract){ # iE <- 1
        iType <- as.character(dt.extract[iE,.SD$type])
        iStart <- dt.extract[iE,.SD$index.start]
        iEnd <- dt.extract[iE,.SD$index.stop]
        if(iType == "chunk"){
            test.export.none <- length(grep(":exports none",file.line[iStart],fixed = TRUE))
            if(rm.export.none && test.export.none==1){
                next
            }
            iStart <- iStart + 1
            iEnd <- iEnd - 1
        }
        if(dt.extract[iE,.SD$export]==FALSE){
            next
        }
        add.before <- switch(iType,
                             "header"=NULL,
                             "chunk"=if(index.chunk){paste0("## chunk ",dt.extract[iE,.SD$index],"")}else{NULL})

        prefix.add <- switch(iType,
                             "header"="## ",
                             "chunk"=NULL)
        toAdd <- paste0(prefix.add,file.line[iStart:iEnd])
        
        add.after <- switch(iType,
                            "header"=if(iE<n.extract&&dt.extract[iE+1,.SD$type]=="chunk"){""}else{NULL},
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
