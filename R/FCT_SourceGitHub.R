# from http://stackoverflow.com/questions/27297484/r-using-rvest-package-instead-of-xml-package-to-get-links-from-url

read_htmlLinks <- function(x, ...){
  require(rvest)
  require(pipeR) # %>>% will be faster than %>%
  
  content <- xml2:::read_html(x,...)
  links <- content %>>% rvest:::html_nodes("a") %>>% rvest:::html_attr("href")
  return(links)
}

read_githubCommit <- function(repo, keep.author = TRUE, keep.name = TRUE, keep.time = TRUE, n.commit = 5, trace = TRUE){
  
  # as in devtools:::github_remote
  meta <- devtools:::parse_git_repo(repo)
  meta <- devtools:::github_resolve_ref(do.call(devtools:::`%||%`, list(meta$ref,"master")), meta)
  if (is.null(meta$username)) {
      stop("Unknown username. Possible mispecification of repo. Have you included the user name?")
  }
 
  #### get all commits
  pathURL <- file.path("https://github.com",meta$username,meta$repo,"commits")
  links <- read_htmlLinks(pathURL)
  pathURL <- file.path(meta$username,meta$repo,"tree")
  index.commits <- grep(pathURL,links, value = FALSE)
  author.commits <- links[index.commits - 2]
  author.commits <- unlist(lapply(strsplit(author.commits, split = "author="),"[[",2))
  
  commits <- links[index.commits]
  commits <- unlist(lapply(strsplit(commits, split = "tree/"),"[[",2))
  n.commits <- min(length(commits),n.commit)
  commits <- commits[1:n.commits]
  
  #### get additional informations
  if(keep.name || keep.time){ 
    
    name.commits <- vector(length = n.commits)
    time.commits <- vector(length = n.commits)
    
    if(trace){
      cat("Find the name of the ",n.commits," last commits \n", sep = "")
      pb <- utils::txtProgressBar(0, n.commits, style = 3)
      on.exit(close(pb))
    }
    
      for(iterC in 1:n.commits){
      pathURL <- file.path("https://github.com",meta$username,meta$repo,"commit",commits[iterC])
      
      # system.time(
        # commit.content <- RCurl:::getURL(pathURL, httpheader = c(name = "description", name = "go-import"))
        commit.content <- read_html(pathURL)
       
        name.commit <- grep("property=\"og:title\"",commit.content %>% html_nodes("meta"), value = TRUE)
        name.commits[iterC] <- strsplit(gsub("<meta content=\"","",name.commit), split = " · ")[[1]][1]
        time.commit <- grep("property=\"og:updated_time\"",commit.content %>% html_nodes("meta"), value = TRUE)
        time.commit <- strsplit(gsub("<meta content=\"","",time.commit), split = "\" property=")[[1]][1]
        time.commits[iterC] <- as.character(as.Date(as.numeric(time.commit)/(24*3600)))
         # commit.content <- substr(RCurl:::getURL(pathURL),
        #                        nchar[1],nchar[2])
        # XML:::htmlTreeParse(commit.content, asText = TRUE)
        # 
        # name.commits[iterC] <- tail(strsplit(split = "<title>",
        #                                      strsplit(split = "GitHub</title>", commit.content)[[1]][1]
        # )[[1]],1)
      if(trace){utils::setTxtProgressBar(pb,iterC)}
    }
  }
 
  
  #### export
  if(keep.name){attr(commits, "name") <- name.commits}
  if(keep.time){attr(commits, "time") <- time.commits}
  if(keep.author){attr(commits, "author") <- author.commits}
  
  return(commits)
}

