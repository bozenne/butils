#' @title links of a url page
#' 
#' @description Get the links of a url page.
#'
#' @param x the url page. Character.
#' @param ... additional arguments to be passed to \code{read_html}.
#' 
#' @details inspired from http://stackoverflow.com/questions/27297484/r-using-rvest-package-instead-of-xml-package-to-get-links-from-url 
#' 
#' @keywords function url
#' @examples 
#' read_htmlLinks("https://cran.r-project.org/")
#' read_htmlLinks("https://github.com/bozenne/BuyseTest")
#' 
#' @export
read_htmlLinks <- function(x, ...){
  # pipeR:::`%>>%` will be faster than rvest:::`%>%`
  
  content <- xml2:::read_html(x,...)
  links <- content %>>% rvest:::html_nodes("a") %>>% rvest:::html_attr("href")
  return(links)
}

#' @title Existing commits
#' 
#' @description Get the existing commits corresponding to a repo on Github
#'
#' @param repo Repository address in the format username/repo
#' @param keep.author should the author be returned
#' @param keep.name should the name of each commit be returned
#' @param keep.time should the time at which the commit was done be returned
#' @param seq.commit the sequence of commit to return
#' @param rev.commit shouidl the order of the commit be reverested (last becomes first).
#' @param trace should the progression in geting the name/time of the commit be displayed
#' @param ... additional arguments to be passed to \code{read_html}.
#' 
#' @details 
#' 
#' @keywords function github
#' @examples 
#' read_githubCommit("bozenne/BuyseTest")
#' read_githubCommit("bozenne/BuyseTest",rev.commit = TRUE)
#' 
#' ## install the first version of the package
#' # devtools:::install_github("bozenne/BuyseTest", ref = "70e5ef6c3c66d44f9227297225e2b9d2fe6acc06")
#' 
#' @export
read_githubCommit <- function(repo, seq.commit = 1:5, rev.commit = FALSE,
                              keep.author = TRUE, keep.name = TRUE, keep.time = TRUE,  trace = TRUE, ...){
  
  # as in devtools:::github_remote
  meta <- devtools:::parse_git_repo(repo)
  meta <- devtools:::github_resolve_ref(do.call(devtools:::`%||%`, list(meta$ref,"master")), meta)
  if (is.null(meta$username)) {
      stop("Unknown username. Possible mispecification of repo. Have you included the user name?")
  }
 
  #### get all commits
  pathURL <- file.path("https://github.com",meta$username,meta$repo,"commits")
  links <- read_htmlLinks(pathURL, ...)
  pathURL <- file.path(meta$username,meta$repo,"tree")
  index.commits <- grep(pathURL,links, value = FALSE)
  author.commits <- links[index.commits - 2]
  author.commits <- unlist(lapply(strsplit(author.commits, split = "author="),"[[",2))
  
  commits <- links[index.commits]
  commits <- unlist(lapply(strsplit(commits, split = "tree/"),"[[",2))
  
  if(rev.commit){commits <- rev(commits)}
  seq.commit <- intersect(seq.commit, seq_len(length(commits)))
  author.commits <- author.commits[seq.commit]
  commits <- commits[seq.commit]
  n.commits <- length(commits)
  
  #### get additional informations
  if(keep.name || keep.time){ 
    
    name.commits <- vector(length = n.commits)
    time.commits <- vector(length = n.commits)
    
    if(trace){
      cat("Find the name of the ",paste(seq.commit, collapse = " "),if(rev.commit){" first"}else{" last"}," commits \n", sep = "")
      pb <- utils::txtProgressBar(0, n.commits, style = 3)
      on.exit(close(pb))
    }
    
      for(iterC in 1:n.commits){
      pathURL <- file.path("https://github.com",meta$username,meta$repo,"commit",commits[iterC])
      
       commit.content <- xml2:::read_html(pathURL, ...)
       
        name.commit <- grep("property=\"og:title\"",commit.content %>>% rvest:::html_nodes("meta"), value = TRUE)
        name.commits[iterC] <- strsplit(gsub("<meta content=\"","",name.commit), split = " · ")[[1]][1]
        time.commit <- grep("property=\"og:updated_time\"",commit.content %>>% rvest:::html_nodes("meta"), value = TRUE)
        time.commit <- strsplit(gsub("<meta content=\"","",time.commit), split = "\" property=")[[1]][1]
        time.commits[iterC] <- as.character(zoo::as.Date(as.numeric(time.commit)/(24*3600))) # must be zoo and not base
      
      if(trace){utils::setTxtProgressBar(pb,iterC)}
    }
  }
 
  
  #### export
  if(keep.name){commits <- cbind(commits, name = name.commits)}
  if(keep.time){commits <- cbind(commits, time = time.commits)}
  if(keep.author){commits <- cbind(commits, author = author.commits)}
  
  return(commits)
}

#' @title Install commit
#' 
#' @description Install a specific commit of a package, possibly under a temporary name
#'
#' @param repo Repository address in the format username/repo
#' @param ref Desired git reference. See the documentation of \code{devtools::install_github}.
#' @param subdir Subdirectory within repo that contains the R package.  See the documentation of \code{devtools::install_github}.
#' @param host GitHub API host to use. See the documentation of \code{devtools::install_github}.
#' @param temporary Should the package be installed under a newname, i.e. adding \code{Tempo} to the name of the package
#' @param force Not sure what it is !!!
#' @param quiet Not sure what it is !!!
#' @param ... additional arguments to be passed to \code{devtools::install_packages} or \code{devtools:::install}
#' 
#' @details 
#' 
#' @keywords function github package
#' @examples 
#' install_githubCommit("bozenne/BuyseTest", ref = "70e5ef6c3c66d44f9227297225e2b9d2fe6acc06")
#' 
#' 
#' @export
install_githubCommit <- function (repo, ref = "master", subdir = NULL, host = "https://api.github.com",
                                  temporary = TRUE, force = FALSE, quiet = FALSE, ...) 
{
  remote <- lapply(repo, devtools:::github_remote, ref = ref, host = host,
                   username = NULL, subdir = subdir)[[1]]
  
  stopifnot(devtools:::is.remote(remote))
  remote_sha <- devtools:::remote_sha(remote)
  package_name <- devtools:::remote_package_name(remote)
  local_sha <- devtools:::local_sha(package_name)
  if (!isTRUE(force) && !devtools:::different_sha(remote_sha = remote_sha, 
                                       local_sha = local_sha)) {
    if (!quiet) {
      message("Skipping install of '", package_name, "' from a ", 
              sub("_remote", "", class(remote)[1L]), " remote,", 
              " the SHA1 (", substr(remote_sha, 1L, 8L), ") has not changed since last install.\n", 
              "  Use `force = TRUE` to force installation")
    }
    return(invisible(FALSE))
  }
  if (devtools:::is_windows && inherits(remote, "cran_remote")) {
    devtools::install_packages(package_name, repos = remote$repos, 
                     type = remote$pkg_type, ..., quiet = quiet)
    return(invisible(TRUE))
  }
  bundle <- devtools:::remote_download(remote, quiet = quiet)
  on.exit(unlink(bundle), add = TRUE)
  source <- devtools:::source_pkg(bundle, subdir = remote$subdir)
  on.exit(unlink(source, recursive = TRUE), add = TRUE)
  metadata <- devtools:::remote_metadata(remote, bundle, source)
  list.files(source)
  
  if(temporary){
    cat("* Installation as a temporary package \n")
    
    # description
    file.description <- readLines(file.path(source,"DESCRIPTION"))
    indexLine <- grep("Package:",file.description)
    packageName <- gsub("Package: ", replacement = "", x = file.description[indexLine])
    file.description[1] <- paste0("Package: ", packageName,"Tempo")
    
    con <- file(file.path(source,"DESCRIPTION")) 
    writeLines(text = file.description, con = con) 
    close(con) 
    
    # namespace
    file.namespace <- readLines(file.path(source,"NAMESPACE"))
    indexLine <- grep("useDynLib",file.namespace)
    if(length(indexLine)!=0){ # cpp links to R
      file.namespace[indexLine] <-  gsub(pattern = paste0("useDynLib(",packageName,")"), 
                                           replacement = paste0("useDynLib(",packageName,"Tempo)"),
                                           x = file.namespace[indexLine], fixed = TRUE)
      
      con <- file(file.path(source,"NAMESPACE")) 
      writeLines(text = file.namespace, con = con) 
      close(con) 
      
      Rcpp::compileAttributes(source)
    }
    
  }

  devtools:::install(source, ..., quiet = quiet, metadata = metadata)
}
