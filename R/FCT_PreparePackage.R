### exportPattern("^[^\\.]") : exporte toutes les fonctions ne commencant pas par un point
### exportPattern("^[[:alpha:]]+") : exporte toutes les fonctions commencant par une lettre 
### "The "Writing R Extensions" manual suggests "the directive exportPattern("^[^\\.]") exports all variables that do not start with a period."

#' @title Integrative function for building packages
#' @description  Convenient way to compile, test or install a package
#' 
#' @param package the name of the package
#' @param version the version of the package
#' @param path.package the position of the directory containing the package
#' @param compileAttributes should the function \code{compileAttributes} be used to enerates the bindings required to call C++ functions from R for functions adorned with the Rcpp::export attribute.
#' @param clearDir should the temporary directory used to unzip the package be removed
#' @param roxygenise should the documentation be generated using \code{roxygenise}
#' @param build should the package be build
#' @param options.build additional options used to build the package
#' @param untar should the package be unzipped
#' @param check should the CRAN test be applied to the package
#' @param options.check additional options used to check the package
#' @param install should the package be installed on the computer.
#' @param options.install additional options used to install the package
#' @param trace How the execution of the function should be traced. If 2, also display the execution of compileAttributes.
#' @param clearAfter Should the directory where the package has been unzipped or tested be removed at the end of the execution.
#' 
#' 
#' @keywords function package
#' 
#' @export
buildPackage <- function(package, version = NULL, path.package = NULL, 
                         compileAttributes = TRUE, clearDir = TRUE, roxygenise = TRUE,
                         build = TRUE, options.build = NULL, 
                         untar = TRUE, 
                         check = FALSE, options.check = NULL, 
                         install = TRUE, options.install = "--build", 
                         trace = 2, clearAfter = TRUE){
  require(Rcpp)
  require(RcppArmadillo)
  require(roxygen2)
  
  if(!is.null(path.package)){
    oldWD <- getwd()
    setwd(path.package)
    on.exit(setwd(oldWD))
  }
  
  if(is.null(version)){
    fieldVersion <- grep("Version", readLines(file.path(path.package, package,"DESCRIPTION")), value = TRUE)
    version <- gsub("Version: ", replacement = "", x = fieldVersion)
  }
  
  packageVersion <- paste(package, version, sep = "_")
  path.Wpackage <- package
  path.WpackageVersion <- packageVersion
  
  if(compileAttributes){
    if(trace)cat(">> compileAttributes \n")
    Rcpp::compileAttributes(pkgdir = path.Wpackage, verbose = (trace>=2))
  } 
  # crlf2lf("testCpp")
  # pdf2qpdf("MRIaggr")
  
  if(clearDir){ 
    if(trace)cat(">> clearDir \n")
    cleanDir(path.WpackageVersion)
    if(trace)cat("\n")
  }
  if(roxygenise){
    if(trace)cat(">> roxygenise \n")
    roxygen2::roxygenise(package)
    if(trace)cat("\n")
  } 
  
  if(build){
    if(trace)cat(">> build \n")
    system(paste0("R CMD build ",options.build," ",package))
    if(trace)cat("\n")
  } 
  
  if(untar){
    if(trace)cat(">> untar \n")
    untar(paste0(packageVersion,".tar.gz"), exdir = paste0(packageVersion))  
    if(trace)cat("\n")
  }
  
  if(check){
    if(trace)cat(">> CRAN check \n")
    system(paste0("R CMD check ",options.check," ",packageVersion,"/",package))
    if(clearAfter){cleanDir(paste0(path.Wpackage,".Rcheck"), test = FALSE)}
    if(trace)cat("\n")
  }
  
  if(install){
    if(trace)cat("* install package \n")
    system(paste0("R CMD INSTALL ",options.check," ",packageVersion,"/",package))
    if(clearAfter){cleanDir(path.WpackageVersion, test = FALSE)}
    if(trace)cat("\n")
  }
  
  
}


#' @title Convert .cpp files with CRLF line ending in LF line ending
#'
#' @references  http://lists.r-forge.r-project.org/pipermail/rcpp-devel/2012-April/003731.html Dirk Eddelbuettel edd at debian.org (Mon Apr 23 18:53:37 CEST 2012)
crlf2lf <- function(dirname){ 
  cat("Convert file ending from CRLF to LF (dir : ",dirname,") \n")
  
  if("src" %in% list.files(dirname) == FALSE){
    cat("\'src\' directory not found \n")
    return(invisible(NULL))
  }
  
  file_cpp <- tools::list_files_with_exts(paste(dirname,"/src",sep=""), exts="cpp", full.names = FALSE)
  
  if(length(file_cpp)==0){
    cat("no \'.cpp\' files founded in the \'src\' directory \n")
    return(invisible(NULL))
  }
  
  n.files <- length(file_cpp)
  cat(n.files,"files founded \n")
  test <- 0
  
  for(iter_file in file_cpp){
    if(test==0){
      valid <- FALSE
      while(valid==FALSE){
        test <- readline(paste("should the file \'",iter_file,"\' be converted ? (1/0) : ",sep=""))
        test <- as.numeric(test)
        
        if(test %% 1 != 0 || test>n.files){
          test <- cat("incorrect answer. Answer should be 1 (yes), 0 (no) or -1 (exit) \n")
        }else{valid <- TRUE}
      }
    }
    
    if(test>=1){
      con <- file(paste(dirname,"/src/",iter_file,sep=""), "rb") 
      bin <- readBin(con, raw(), 100000) 
      bin <- bin[ which(bin != "0d") ] 
      close(con) 
      
      Sys.sleep(1) 
      
      con <- file(paste(dirname,"/src/",iter_file,sep=""), "wb") 
      writeBin(bin, con) 
      close(con) 
      cat("done \n") 
      
      test <- test-1
    }
    if(test<0){
      cat("exit \n") 
      break}
  }
}

#' @title Remove all files from a given directory
#'
#' @param dirname the path leading to the directory
#' @param test should the user be asked whether he really wants to remove the directory
cleanDir <- function(dirname,test=TRUE){
  
  dirname_split <- strsplit(dirname,split="/")[[1]]
  if(length(dirname_split)>1){
    path <- paste(dirname_split[-length(dirname_split)],collapse="/")
    dirname_1 <- dirname_split[length(dirname_split)]
  }else{
    path <- "."
    dirname_1 <- dirname_split
  }
  
  list_files <- list.files(dirname,recursive = TRUE)
  n.files <- length(list_files)
  if(dirname_1 %in% list.files(path)){
    if(test==TRUE){
      test <- readline(paste("should the directory \'",dirname,"\' be removed ? (1/0) : \n",
                             "(",n.files," files )",sep=""))
    }else{
      test <- 1
    }
    
    if(test==1){
      unlink(dirname,recursive=TRUE)
      cat("done \n") 
    }
  }else{cat("directory not in the current working directory \n")}
  
  return(invisible(list_files))
}

#' @title Convert from pdf to qpdf
#'
pdf2qpdf <- function(dirname,path_qpdf="../qpdf-5.1.2-bin-mingw32/qpdf-5.1.2/bin/qpdf"){
  cat("Convert pdf to qpdf (dir : ",dirname,") \n")
  
  if("doc" %in% list.files(dirname) == FALSE){
    cat("\'doc\' directory not found \n")
    return(invisible(NULL))
  }
  
  file_pdf <- list.files(paste(dirname,"/doc",sep=""))
  file_pdf <- file_pdf[grep(file_pdf,pattern=".pdf",fixed=TRUE)]
  
  if(length(file_pdf)==0){
    cat("no \'.pdf\' files founded in the \'doc\' directory \n")
    return(invisible(NULL))
  }
  
  n.files <- length(file_pdf)
  cat(n.files,"files founded \n")
  test <- 0
  
  for(iter_file in file_pdf){
    if(test<=1){
      test <- readline(paste("should the file \'",iter_file,"\' be converted ? (1/0) : ",sep=""))
    }
    if(test>=1){
      system(paste(path_qpdf," ",dirname,"/doc/",iter_file," ",dirname,"/doc/",iter_file,sep=""))
      cat("done \n") 
    }
    if(test<0){
      cat("exit \n") 
      break}
  }
}

#' @title Print all the arguments from the Rd file
#'
printLsArgs <- function(dir,write.file=FALSE,filename="List_of_arguments",trace=TRUE,path=NULL){
  
  file_Rd <- list.files(dir)
  file_Rd <- file_Rd[grep(file_Rd,pattern=".Rd",fixed=TRUE)]
  
  if(length(file_Rd)==0){
    cat("no \'.Rd\' files founded in the \'",dir,"\' directory \n",sep="")
    return(invisible(NULL))
  }
  
  n.files <- length(file_Rd)
  if(trace==TRUE){cat(n.files,"files founded \n")}
  
  Margs.files <- matrix(NA,nrow=n.files,ncol=3)
  rownames(Margs.files) <- file_Rd
  colnames(Margs.files) <- c("n.args","line.begin","line.end")
  
  ls.explain <- list()
  
  
  doubleBar <- gsub(pattern=" ",replacement="",x="\\ ")
  singleBar <- gsub(pattern=" ",replacement="",x="\ ")
  
  for(iter_file in 1:n.files){
    
    # load file
    file.size <-   file.info(paste(dir,file_Rd[iter_file],sep=""))$size
    con <- file(paste(dir,file_Rd[iter_file],sep=""), "rb") 
    file.character <- readChar(con, file.size*10) 
    close(con)
    
    # split file in lines
    file.lines <- unlist(strsplit(file.character,split="\n",fixed=TRUE))
    
    # rm \r
    test.slashR <- length(grep(pattern="\r",x=file.lines))>0
    file.lines <- gsub(pattern="\r",replacement="",file.lines)
    
    #     browser(expr=iter_file==62)
    # find argument section
    file.lines.spaceRemoved <- gsub(" ","",file.lines)
    line.begin <- which(file.lines.spaceRemoved=="\\arguments{")+1
    if(length(line.begin)==0){ # if no argument section
      if(trace){cat(iter_file,"- file \"",file_Rd[iter_file],"\" ",if(test.slashR){"(rm \\r)"}," : no arguments section \n",sep="")}
      next
    }
    line.end <- setdiff(which(file.lines.spaceRemoved=="}"),1:line.begin)[1]-1
    
    Margs.files[file_Rd[iter_file],"line.begin"] <- line.begin
    Margs.files[file_Rd[iter_file],"line.end"] <- line.end
    
    
    # find non empty lines of the argument section
    arg.lines <- (line.begin:line.end)[nchar(file.lines.spaceRemoved[line.begin:line.end])>0]
    
    # find extra lines
    Rarg.lines <- arg.lines[grep(pattern="\\item{",x=file.lines[arg.lines],fixed=TRUE)]
    extra.lines <- setdiff(arg.lines,Rarg.lines)
    Margs.files[file_Rd[iter_file],"n.args"] <- n.args <- length(Rarg.lines)  
    
    ls.extra.lines <- sapply(1:n.args,function(x){
      extra.lines[(extra.lines > Rarg.lines[x])*(extra.lines < c(Rarg.lines,Inf)[x+1])==1]
    })
    
    # display
    if(trace){cat(iter_file,"- file \"",file_Rd[iter_file],"\" ",if(test.slashR){"(rm \\r)"}," : ",n.args," arguments (line ",line.begin," to ",line.end,")\n",sep="")}
    
    
    for(iter_arg in 1:n.args){
      
      iter_line <- Rarg.lines[iter_arg]
      
      if(length(ls.extra.lines[[iter_arg]])>0){ # merge extra lines
        file.lines[iter_line] <- paste(file.lines[iter_line],
                                       file.lines[ls.extra.lines[[iter_arg]]],sep="\n") 
      }
      
      # select argument and explaination
      arg_tempo <- strsplit(strsplit(file.lines[iter_line],split="\\item{",fixed=TRUE)[[1]][2],split="}{",fixed=TRUE)[[1]][1]
      explain_tempo <- paste("Line ",iter_line," : ",file.lines[iter_line]," # ",file_Rd[iter_file],"\n",sep="")
      
      if(arg_tempo %in% colnames(Margs.files)==FALSE){ # add new argument
        Margs.files <- cbind(Margs.files,NA)
        colnames(Margs.files)[ncol(Margs.files)] <- arg_tempo
        ls.explain[[arg_tempo]] <- explain_tempo
      }else{
        ls.explain[[arg_tempo]] <- c(ls.explain[[arg_tempo]],explain_tempo)
      }
      
      Margs.files[file_Rd[iter_file],arg_tempo] <- T
      
    }
    
    
  }
  
  # rank by alphabetic order
  ls.explain <- ls.explain[sort(names(ls.explain))]
  
  
  # write file
  if(write.file==TRUE){
    
    if(trace==TRUE){cat("creating and writting \'",filename,"\' file in ",if(!is.null(path)){path}else{"the current directory"}," \n",sep="")}
    
    ls.text <- sapply(1:length(ls.explain),function(x){
      paste("#### ",x,"- ",names(ls.explain)[x]," #### \n",paste(ls.explain[[x]],collapse=""),"\n",sep="")
    })
    
    text <- c(as.character(unlist(ls.text)),
              "\n \n################################# ",n.files," files .Rd #################################\n",
              paste(sort(rownames(Margs.files)),"\n",sep="")
    )
    
    file.create <- file(paste(path,filename,sep=""), "wb")
    writeLines(text=text,con=file.create, sep="")
    close(file.create)
  }
  
  # export
  return(list(Margs.files=Margs.files,
              ls.explain=ls.explain))
  
}



