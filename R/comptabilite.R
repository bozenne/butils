### comptabilite.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec  2 2017 (12:29) 
## Version: 
## Last-Updated: feb  6 2018 (17:15) 
##           By: Brice Ozenne
##     Update #: 390
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
#' @param involved [vector of characters] who was involved in the activity?
#' @param type [character] the name of the activity.
#' @param date [date] the date at which the activity happen.
#' @param note [character] additional text.
#' @param value [numeric] a named vector describing name paid what.
#' @param label [label] a label associated with the activity. When
#' displayed activities are sorted by labels.
#' @param ... ignored.
#' 
#' @examples
#' #### Create an account ####
#' myAcc <- createAccount()
#'
#' #### Ease input with short names ####
#' addNickName(myAcc) <- list("Brice Ozenne" = "Brice",
#'                            "Sebastian Holst" = c("Seb","Sebastian"))
#'
#' #### Add activity paid by one person for others ####
#' addActivity(myAcc, involved = c("Brice","Sebastian"), type = "Court") <- c("Sebastian" = 200)
#'
#' #### Add expense, i.e. activity paid by nobody ####
#' addActivity(myAcc, involved = c("Brice","Sebastian"), type = "Shuttlecock") <- 100
#'
#' #### Add contribution from one member ####
#' addActivity(myAcc, type = "Shuttlecock") <- c("Brice" = 100)
#'
#' ####
#' myAcc
#' @export
`addActivity<-` <-
    function(object, ..., value) UseMethod("addActivity<-")

#' @rdname addActivity
#' @export
"addActivity<-.butilsAccount" <- function(object,
                                          involved = NULL,
                                          type = as.character(NA),
                                          date = as.Date(NA),
                                          note = "",
                                          label = NULL,
                                          ...,
                                          value){

    ### ** normalize arguments
    if(identical(names(value),"")){
        names(value) <- NULL
    }
    who.paid <- names(value)
    
    test.cost <- !is.null(involved)    
    test.paid <- !is.null(who.paid)
    
    ### ** check consistency of the arguments
    if(length(value)!=1){
        stop("Argument \'value\' must have length 1 \n")
    }
    
    if(is.null(involved) && is.null(who.paid)){
        stop("argument \'involved\' must be specified or argument \'value\' must be named \n")
    }
    
    if(!is.null(involved) && !is.character(involved)){
        stop("argument \'involved\' should be a character vector \n")
    }
    
    if(!is.character(type)){
        stop("argument \'type\' should be a character vector \n")
    }
    
    if("Date" %in% class(date) == FALSE){
        stop("argument \'Date\' should inherit from the class \"Date\" \n")
    }
    
    ### ** Convert nick names to real names
    if(!is.null(object$nickName)){

        for(iName in names(object$nickName)){ ## iName <- names(object$nickName)[1]
            nickNames <- paste0("^",paste(object$nickName[[iName]], collapse = "$|^"),"$")

            if(!is.null(involved)){
                index <- grep(pattern = nickNames, x = involved)
                if(length(index)>0){
                    involved[index] <- iName
                }
            }

            if(!is.null(who.paid)){
                index <- grep(pattern = nickNames, x = who.paid)
                if(length(index)>0){
                    who.paid <- iName
                }
            }
        }

    }

    ### ** prepare
    
    ### *** payment
    all.involved <- union(involved, who.paid)
    vec.paid <- rep(0, length(all.involved))

    if(test.paid){
        ## Contribution of one of the member
        vec.paid[who.paid == all.involved] <- as.double(value)
    }else{
        ## No contribution from the members
        ## e.g. 10 shuttlecocks have been used
        ## vec.paid stay at 0 for everybody
    }

    
    
    ### *** cost
    if(test.cost){
        ## Cost related to the activity
        ## the total cost will be shared among the people involved (see below)
        cost <- as.double(value)
    }else{
        ## Contribution of one of the member (nobody is involved)
        ## e.g. Brice buy shuttlecock for 100 kr. They have not yet been used thought.
        ##      It does not (yet) corresponds to an expanse for anybody
        ##      This is why total cost is set to 0
        cost <- 0
    }
    
    ### *** Unique identifier for each entry
    if(is.null(object$table)){
        id.entry <- 1
    }else{
        id.entry <- max(object$table$id.entry)+1
    }

    if(is.null(label)){
        label <- id.entry
    }

    ### ** fill the table
    newtable <- data.table(paid = vec.paid,
                           name = all.involved,
                           date = date,
                           note = note,
                           type = type,
                           label = label,
                           id.entry = id.entry)

    newtable[, c("total.cost","n.participant","participant.cost") := 0]
    newtable[name %in% involved, "total.cost" := cost]
    newtable[name %in% involved, "n.participant" := .N]
    newtable[name %in% involved, "participant.cost" := .SD$total.cost/.SD$n.participant]

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
#' @param print [logical] should the summary be printed in the console.
#' @param detail [logical] should the balance per individual (1) and by activity (2) be displayed.
#' @param keep.cols [character vector] the columns to be displayed in the detail.
#' @param digit [integer] the number of decimal to be displayed.
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
                                  detail = 0:2,
                                  keep.cols = c("paid","date","type","total.cost","participant.cost"),
                                  digit = 1,
                                  ...){

    balance <- label <- paid <- spent <- NULL ## [:CRAN checks] data.table
    
    ### ** Count
    if(!is.null(object$table)){
        text.cat <- "#### balance ####\n"
        
        balance.print <- object$table[,list(paid = sum(.SD$paid), spent = sum(.SD$participant.cost)),
                                      by = "name",
                                      .SDcols = c("paid","participant.cost")]
        balance.print[, c("balance") :=  .SD$paid - .SD$spent,
                      .SDcols = c("paid","spent")]

        setkeyv(object$table, c("date","id.entry"))
        tempo1 <- object$table[,list(list(.SD)), .SDcols = keep.cols, by = "name"]
        detail1.print <- setNames(tempo1[[2]],tempo1[[1]])

        setkeyv(object$table, c("label", "paid"))

        Ulab <- as.character(object$table[,1,by="label"][,label])
        n.label <- length(Ulab)
        detail2.print <- vector(mode = "list", length = n.label)
        names(detail2.print) <- Ulab
        for(iLab in 1:n.label){ ## iLab <- 1  
            detail2.print[[Ulab[iLab]]] <- object$table[label == Ulab[iLab],.SD[1], .SDcols = setdiff(keep.cols, "paid"), by = "id.entry"]                   
        }
    }else{
        text.cat <- "the account is empty \n"
        detail <- NULL
        balance.print <- NULL
        detail1.print <- NULL
        detail2.print <- NULL       
    }
    
    ### ** Display

    
    if(!is.null(detail)){
        if(0 %in% detail){
            cat(text.cat)
            balance.print[, paid := round(paid, digits = digit)]
            balance.print[, spent := round(spent, digits = digit)]
            balance.print[, balance := round(balance, digits = digit)]
            balance.print <- as.data.frame(balance.print)
            print(balance.print)
        }
        
        if(1 %in% detail){
            if(0 %in% detail){
                cat("\n\n")
            }
            
            cat("#### detail of the spending by individual ####")
            detail1.print <- lapply(detail1.print, as.data.frame)
            for(iName in names(detail1.print)){
                cat("\n ")
                iDF <- detail1.print[[iName]]
                cat("*** name: ",iName,"\n", sep = "")
                iDF$paid <- round(iDF$paid, digits = digit)
                iDF$total.cost <- round(iDF$total.cost, digits = digit)
                iDF$participant.cost <- round(iDF$participant.cost, digits = digit)
                print(iDF)
            }
        }

        if(2 %in% detail){
            if(any(0:1 %in% detail)){
                cat("\n\n")
            }
            
            cat("#### detail of the spending by activity ####")
            detail2.print <- lapply(detail2.print, as.data.frame)
            for(iAct in names(detail2.print)){
                cat("\n ")
                iDF <- detail2.print[[iAct]]
                cat("*** activity: ",iAct,"\n", sep = "")
                iDF$total.cost <- round(iDF$total.cost, digits = digit)
                iDF$participant.cost <- round(iDF$participant.cost, digits = digit)
                print(iDF)
            }
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
