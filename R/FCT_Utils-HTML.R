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
