### comptabilite.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec  2 2017 (12:29) 
## Version: 
## Last-Updated: dec 20 2017 (18:36) 
##           By: Brice Ozenne
##     Update #: 154
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


## * createAccount
#' @title Create an empty account
#' @description Create an empty account.
#' 
#' @export
createAccount <- function(){
    object <- list(table = NULL, nickName = NULL)
    class(object) <- "butilsAccount"
    return(object)
}
## * addNickName
#' @title Add nick name to the account.
#' @description Add nick name to the account.
#' @name addNickName
#' 
#' @param object the account.
#' @param value a named vector.
#' 
#' @examples
#' myAcc <- createAccount()
#' addNickName(myAcc) <- list("Brice Ozenne" = "Brice",
#'                            "Sebastian Holst" = c("Seb","Sebastian"))
#' myAcc
#' @export
`addNickName<-` <-
    function(object, value) UseMethod("addNickName<-")

#' @rdname addNickName
#' @export
"addNickName<-.butilsAccount" <- function(object,
                                          value){

### ** Update
    if(is.null(object$nickName)){
        object$nickName <- value    
    }else{
        object$nickName <- c(object$nickName,value)
    }

### ** Check consistency
    vec.name1 <- names(object$nickName)
    if(any(duplicated(vec.name1))){
        stop("duplicated full names \n",
             "\"",paste(unique(vec.name1[duplicated(vec.name1)]), collapse = "\" \""),"\"\n")
    }
    vec.name2 <- as.character(unlist(object$nickName))
    if(any(duplicated(vec.name2))){
        stop("duplicated nick names \n",
             "\"",paste(unique(vec.name1[duplicated(vec.name2)]), collapse = "\" \""),"\"\n")
    }

    if(any(vec.name1 %in% vec.name2)){
        stop("some full names and some nick names coincide \n",
             "\"",paste(unique(vec.name1[vec.name1 %in% vec.name2]), collapse = "\" \""),"\"\n")
    }

### ** export    
    return(object)
}

## * addActivity
#' @title Add a new entry to the account.
#' @description Add a new entry to the account.
#' @name addActivity
#' 
#' @param object the account.
#' @param involved who was involved in the activity?
#' @param type a character string describing the activity.
#' @param date the date at which the activity happen.
#' @param value a named vector describing name paid what.
#' @param ... ignored.
#' 
#' @examples
#' myAcc <- createAccount()
#' addNickName(myAcc) <- list("Brice Ozenne" = "Brice",
#'                            "Sebastian Holst" = c("Seb","Sebastian"))
#' addActivity(myAcc, involved = c("Brice","Sebastian"), type = "Court") <- c("Sebastian" = 95)
#' myAcc
#' @export
`addActivity<-` <-
    function(object, ..., value) UseMethod("addActivity<-")

#' @rdname addActivity
#' @export
"addActivity<-.butilsAccount" <- function(object,
                                          involved,
                                          type = as.character(NA),
                                          date = as.Date(NA),
                                          note = "",
                                          ...,
                                          value){

### ** check consistency of the arguments
    if(!is.character(involved)){
        stop("argument \'involved\' should be a character vector \n")
    }    
    if(!is.character(type)){
        stop("argument \'type\' should be a character vector \n")
    }
    if("Date" %in% class(date) == FALSE){
        stop("argument \'Date\' should inherit from the class \"Date\" \n")
    }
    if(is.null(names(value))){
        stop("argument \'value\' should be named to know name has paid what \n")
    }
    if(any(names(value) %in% involved == FALSE)){
        stop("names of argument \'value\' are not consistent with argument \'involved\' \n")
    }

### ** Convert to real name
    name.value <- names(value)
    if(!is.null(object$nickName)){

        for(iName in names(object$nickName)){
            nickNames <- paste(object$nickName[[iName]], collapse = "|")

            index <- grep(pattern = nickNames, x = involved)
            if(length(index)>1){
                stop("Two names matches the same nickName \n",
                     "Full name: \"",iName,"\"\n",
                     "Proposed name: \"",paste(involved[index], collapse = "\" \""),"\" \n")
            }            
            involved[index] <- iName

            index <- grep(pattern = nickNames, x = name.value)
            if(length(index)>1){
                stop("Two names matches the same nickName \n",
                     "Full name: \"",iName,"\"\n",
                     "Proposed name: \"",paste(name.value[index], collapse = "\" \""),"\" \n")
            }            
            name.value[index] <- iName
        }



    }

### ** fill the table
    value.full <- setNames(rep(0, length(involved)), involved)
    value.full[name.value] <- as.double(value)

    if(is.null(object$table)){
        label <- 1
    }else{
        label <- max(object$table$label)+1
    }
    
    
    newtable <- data.table(paid = as.double(value.full),
                           name = involved,
                           date = date,
                           note = note,
                           type = type,
                           label = label)
    newtable[, c("total.price") := sum(.SD$paid), .SDcols = "paid"]
    newtable[, c("n.participant") := .N]
    newtable[, c("participant.price") := .SD$total.price/.SD$n.participant,
             .SDcols = c("total.price","n.participant")]
    
    object$table <- rbind(object$table,
                          newtable)

### ** return
    return(object)
}

## * summary/print
#' @title Summarizing an account.
#' @description Summarizing an account.
#' @name summary.butilsAccount
#' 
#' @param x an object of class butilsAccount.
#' @param object an object of class butilsAccount.
#' @param print should the summary be printed in the console.
#' @param detail should the balance per individual (1) and by activity (2) be displayed.
#' @param keep.cols the columns to be displayed in the detail.
#' @param ... ignored
#' 

#' @rdname summary.butilsAccount
#' @method print butilsAccount
#' @export
print.butilsAccount <- function(x, ...){
    out <- summary(x, print = TRUE, detail = FALSE)
    return(invisible(out))
}

#' @rdname summary.butilsAccount
#' @method summary butilsAccount
#' @export
summary.butilsAccount <- function(object,
                                  print = TRUE,
                                  detail = 1:2,
                                  keep.cols = c("paid","date","type","total.price","participant.price"),
                                  ...){

    ### ** Count
    if(!is.null(object$table)){
        text.cat <- "#### balance ####\n"
        
        balance.print <- object$table[,list(paid = sum(.SD$paid), spent = sum(.SD$participant.price)),
                                      by = "name",
                                      .SDcols = c("paid","participant.price")]
        balance.print[, c("balance") :=  .SD$paid - .SD$spent,
                      .SDcols = c("paid","spent")]
    
    tempo1 <- object$table[,list(list(.SD)), .SDcols = keep.cols, by = "name"]
    detail1.print <- setNames(tempo1[[2]],tempo1[[1]])

    tempo2 <- object$table[,list(list(.SD)), .SDcols = keep.cols, by = "label"]
    detail2.print <- setNames(tempo2[[2]],tempo2[[1]])
    }else{
        text.cat <- "the account is empty \n"
        balance.print <- NULL
        detail <- NULL
        detail1.print <- NULL
        detail2.print <- NULL
    }
    
    ### ** Display

    cat(text.cat)
    if(!is.null(balance.print)){
        print(balance.print)
    }
    
    if(length(detail)>0){
        if(1 %in% detail){
            cat("\n\n")
            cat("#### detail of the spending by individual ####\n")
            print(detail1.print[])
        }

        if(2 %in% detail){
            cat("\n")
            cat("#### detail of the spending by activity ####\n")
            print(detail2.print[])
        }
    }

    ### ** Export
    out <- list(balance = balance.print,
                detail.indiv = detail1.print,
                detail.activity = detail2.print)

    return(invisible(out))

}

##----------------------------------------------------------------------
### comptabilite.R ends here
