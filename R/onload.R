### onload.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: sep 16 2017 (23:42) 
## Version: 
## last-updated: dec  6 2019 (18:32) 
##           By: Brice Ozenne
##     Update #: 12
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## Functions not exported by boot
boot_isMatrix <- get("isMatrix", envir = asNamespace("boot"), inherits = FALSE)
boot_index.array <- get("index.array", envir = asNamespace("boot"), inherits = FALSE)
boot_boot.return <- get("boot.return", envir = asNamespace("boot"), inherits = FALSE)

lavaSearch2_evalInParentEnv <- get("evalInParentEnv", envir = asNamespace("lavaSearch2"), inherits = FALSE)

estimate.lvm <- get("estimate.lvm", envir = asNamespace("lava"), inherits = FALSE) ## needed for bootReg.lvm since estimate.lvm is saved in the call


validCharacter <- get("validCharacter", envir = asNamespace("BuyseTest"), inherits = FALSE) 
validNames <- get("validNames", envir = asNamespace("BuyseTest"), inherits = FALSE) 
validPath <- get("validPath", envir = asNamespace("BuyseTest"), inherits = FALSE) 
#----------------------------------------------------------------------
### onload.R ends here
