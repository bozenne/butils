### sinkDirectory.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 22 2018 (09:24) 
## Version: 
## Last-Updated: mar 22 2018 (09:37) 
##           By: Brice Ozenne
##     Update #: 6
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

##' @title Import All .rds Files in a Directory.
##' @description Import all .rds files in a directory.
##'
##' @param path [character] path to the directory.
##' @param string.keep [regexp] character string indicating files to import. 
##' @param string.exclude [regexp] character string indicating files to ignore.  
##' @param addMissingCol  [logical] if a dataset does not have the same columns as
##' the other, the necessary empty columns are added to it with \code{NA} as values.
##' @param trace [logical] should a progress bar be displayed tracking how many files
##' have been imported.
##'
##' @details The function first read the name of all the files in the directory.
##' Then if the argument \code{string.keep} is not \code{NULL}, it only retains the
##' files whose name contain \code{string.keep}.
##' Then if the argument \code{string.exclude} is not \code{NULL}, it excludes the
##' files whose name contain \code{string.exclude}.
##'
##' Each file must contain a \code{data.table} object with the same columns,
##' so that they can be combined.
##' @return A \code{data.table} object.
##' @author Brice Ozenne
##' @export
sinkDirectory <- function(path, string.keep = NULL, string.exclude = NULL,
                          addMissingCol = FALSE,
                          trace = TRUE){
    allFiles <- list.files(path)
    index.file <- 1:length(allFiles)

### ** subset files
    if(!is.null(string.keep)){
        index.file <- intersect(index.file,grep(string.keep,allFiles,fixed = TRUE))
    }
    if(!is.null(string.exclude)){
        index.file <- setdiff(index.file,grep(string.exclude,allFiles,fixed = TRUE))
    }

    n.files <- length(index.file)

    #### ** merge files
    dt.merge <- NULL
    if(trace){
        cat("read ",n.files," files \n", sep = "")
        pb <- txtProgressBar(max = n.files)
    }
    for(iFile in 1:n.files){ # iFile <- 1

        ### *** read file 
        iFilename <- allFiles[index.file[iFile]]
        iRes <- readRDS(file = file.path(path,iFilename))
        iRes[,iFile := iFilename]

        ### *** solve pb with missing columns        
        if(!is.null(dt.merge) && addMissingCol==TRUE && NCOL(dt.merge)!=NCOL(iRes)){

            missing <- setdiff(names(dt.merge),names(iRes))
            if(!is.null(missing)){
                for(iMiss in missing){ ## iMiss <- missing[1]
                    vec.tempo <- dt.merge[1,.SD,.SDcols = iMiss][[1]]
                    vec.tempo[1] <- NA                    
                    iRes[, c(iMiss) := vec.tempo[1]]                    
                }
            }
            missing <- setdiff(names(iRes),names(dt.merge))
            if(!is.null(missing)){                
                for(iMiss in missing){ ## iMiss <- missing[1]
                    vec.tempo <- iRes[1,.SD,.SDcols = iMiss][[1]]
                    vec.tempo[1] <- NA                    
                    dt.merge[, c(iMiss) := vec.tempo[1]]
                }
            }

        }

        ### *** merge
        dt.merge <- rbind(dt.merge,iRes)
        if(trace){setTxtProgressBar(pb,iFile)}
    }
    if(trace){close(pb)}

#### ** export
    return(dt.merge)
}

##----------------------------------------------------------------------
### sinkDirectory.R ends here
