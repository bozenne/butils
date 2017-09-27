### boot2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: sep 16 2017 (23:37) 
## Version: 
## last-updated: sep 16 2017 (23:50) 
##           By: Brice Ozenne
##     Update #: 17
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * bootPb - Documentation
##' @title Copy of \code{boot::boot} with progress bar
##' @name bootPb
##' 
##' @param data see \code{boot::boot}
##' @param statistic see \code{boot::boot}
##' @param R see \code{boot::boot}
##' @param sim see \code{boot::boot}
##' @param stype see \code{boot::boot}
##' @param strata see \code{boot::boot}
##' @param L see \code{boot::boot}
##' @param m see \code{boot::boot}
##' @param weights see \code{boot::boot}
##' @param ran.gen see \code{boot::boot}
##' @param mle see \code{boot::boot}
##' @param simple see \code{boot::boot}
##' @param ... see \code{boot::boot}
##' @param parallel see \code{boot::boot}
##' @param ncpus see \code{boot::boot}
##' @param cl see \code{boot::boot}
##'
##' @details Copied from https://stackoverflow.com/questions/37673931/add-a-progress-bar-to-boot-function-in-r
##' 
##' @return see \code{boot::boot}
##'
##' @examples
##' 
##' m.boot <- function(data, indices) {
##'   d <- data[indices]
##'   mean(d, na.rm = T) 
##' }
##' tot_rep <- 200
##' v1 <- rnorm(1000)
##' b <- bootPb(v1, m.boot, R = tot_rep)

## * bootPb - Code
##' @rdname bootPb
##' @export
bootPb <- function (data, statistic, R, sim = "ordinary", stype = c("i", 
                                                                    "f", "w"), strata = rep(1, n), L = NULL, m = 0, weights = NULL, 
                    ran.gen = function(d, p) d, mle = NULL, simple = FALSE, ..., 
                    parallel = c("no", "multicore", "snow"), ncpus = getOption("boot.ncpus", 
                                                                               1L), cl = NULL) 
{
    call <- match.call()
    stype <- match.arg(stype)
    if (missing(parallel)) 
        parallel <- getOption("boot.parallel", "no")
    parallel <- match.arg(parallel)
    have_mc <- have_snow <- FALSE
    if (parallel != "no" && ncpus > 1L) {
        if (parallel == "multicore") 
            have_mc <- .Platform$OS.type != "windows"
        else if (parallel == "snow") 
            have_snow <- TRUE
        if (!have_mc && !have_snow) 
            ncpus <- 1L
        loadNamespace("parallel")
    }
    if (simple && (sim != "ordinary" || stype != "i" || sum(m))) {
        warning("'simple=TRUE' is only valid for 'sim=\"ordinary\", stype=\"i\", n=0', so ignored")
        simple <- FALSE
    }
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
        runif(1)
    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    n <- NROW(data)
    if ((n == 0) || is.null(n)) 
        stop("no data in call to 'boot'")
    temp.str <- strata
    strata <- tapply(seq_len(n), as.numeric(strata))
    t0 <- if (sim != "parametric") {
              if ((sim == "antithetic") && is.null(L)) 
                  L <- empinf(data = data, statistic = statistic, stype = stype, 
                              strata = strata, ...)
              if (sim != "ordinary") 
                  m <- 0
              else if (any(m < 0)) 
                  stop("negative value of 'm' supplied")
              if ((length(m) != 1L) && (length(m) != length(table(strata)))) 
                  stop("length of 'm' incompatible with 'strata'")
              if ((sim == "ordinary") || (sim == "balanced")) {
                  if (boot_isMatrix(weights) && (nrow(weights) != length(R))) 
                      stop("dimensions of 'R' and 'weights' do not match")
              }
              else weights <- NULL
              if (!is.null(weights)) 
                  weights <- t(apply(matrix(weights, n, length(R), 
                                            byrow = TRUE), 2L, normalize, strata))
              if (!simple) 
                  i <- boot_index.array(n, R, sim, strata, m, L, weights)
              original <- if (stype == "f") 
                              rep(1, n)
                          else if (stype == "w") {
                              ns <- tabulate(strata)[strata]
                              1/ns
                          }
                          else seq_len(n)
              t0 <- if (sum(m) > 0L) 
                        statistic(data, original, rep(1, sum(m)), ...)
                    else statistic(data, original, ...)
              rm(original)
              t0
          }
          else statistic(data, ...)
    pred.i <- NULL
    fn <- if (sim == "parametric") {
              ran.gen
              data
              mle
              function(r) {
                  dd <- ran.gen(data, mle)
                  statistic(dd, ...)
              }
          }
          else {
              if (!simple && ncol(i) > n) {
                  pred.i <- as.matrix(i[, (n + 1L):ncol(i)])
                  i <- i[, seq_len(n)]
              }
              if (stype %in% c("f", "w")) {
                  f <- freq.array(i)
                  rm(i)
                  if (stype == "w") 
                      f <- f/ns
                  if (sum(m) == 0L) 
                      function(r) statistic(data, f[r, ], ...)
                  else function(r) statistic(data, f[r, ], pred.i[r, 
                                                                  ], ...)
              }
              else if (sum(m) > 0L) 
                  function(r) statistic(data, i[r, ], pred.i[r, ], 
                                        ...)
              else if (simple) 
                  function(r) statistic(data, boot_index.array(n, 1, sim, 
                                                               strata, m, L, weights), ...)
              else function(r) statistic(data, i[r, ], ...)
          }
    RR <- sum(R)
    res <- if (ncpus > 1L && (have_mc || have_snow)) {
               if (have_mc) {
                   parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus)
               }
               else if (have_snow) {
                   list(...)
                   if (is.null(cl)) {
                       cl <- parallel::makePSOCKcluster(rep("localhost", 
                                                            ncpus))
                       if (RNGkind()[1L] == "L'Ecuyer-CMRG") 
                           parallel::clusterSetRNGStream(cl)
                       res <- parallel::parLapply(cl, seq_len(RR), fn)
                       parallel::stopCluster(cl)
                       res
                   }
                   else parallel::parLapply(cl, seq_len(RR), fn)
               }
           }
           else pbapply::pblapply(seq_len(RR), fn) #### changed !!!
    t.star <- matrix(, RR, length(t0))
    for (r in seq_len(RR)) t.star[r, ] <- res[[r]]
    if (is.null(weights)) 
        weights <- 1/tabulate(strata)[strata]
    boot_boot.return(sim, t0, t.star, temp.str, R, data, statistic, 
                     stype, call, seed, L, m, pred.i, weights, ran.gen, mle)
}


#----------------------------------------------------------------------
### boot2.R ends here
