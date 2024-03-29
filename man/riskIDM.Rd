% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/riskIDM.R
\name{riskIDM}
\alias{riskIDM}
\title{Estimate an Illness Death Model}
\usage{
riskIDM(
  formula,
  data,
  PH,
  time = NULL,
  intervention = NULL,
  var.id,
  var.time,
  var.type,
  start.type = NULL,
  n.boot = 0,
  type.ci = "gaussian",
  level = 0.95,
  cpus = 1,
  seed = NULL,
  keep.indiv = FALSE,
  trace = TRUE
)
}
\arguments{
\item{formula}{[formula] A formula indicating baseline covariates on the right-hand side.}

\item{data}{[data.frame] dataset}

\item{PH}{[logical] should the ratio between the hazard with and without exposure be time independent?}

\item{time}{[numeric vector, >=0] times at which the occupancy probability should be estimated.}

\item{intervention}{[list of matrix] list where each element is a matrix used to deduce the intervention hazards by premultiplying the estimated hazard.}

\item{var.id}{[character] name of the column containing the subject id, i.e. unique identifier for each line.}

\item{var.time}{[character vector of length 2] name of the columns containing the time variables, i.e. time at which each type of event happen (intermediate or absorbing).
If an intermediate event does not occur (e.g. no switch of treatment) then the time variable should be set to the end of follow-up time.}

\item{var.type}{[character vector of length 2] name of the columns containing event type indicator.
The first event type indicator can be categorical (multiple intermediate states) but the last one should be binary.}

\item{start.type}{[character] starting state. Deduced from \code{var.type} if left unspecified.}

\item{n.boot}{[interger, >=0] When strictly positive a non-parametric bootstrap is performed to quantify the uncertainty of the estimates.
The value then indicates the number of bootstrap samples.}

\item{type.ci}{[character] Method used to evaluate confidence intervals from bootstrap samples: under normality assumption (\code{"gaussian"}) or non-parametric (\code{"percentile"}).
The later typically require more bootstrap samples to provide reliable results (in term of Monte Carlo error).}

\item{level}{[numeric, 0.1] Confidence level for the confidence interval.}

\item{cpus}{[integer, >0] the number of CPU to use when doing bootstrap resampling. Default value is 1.}

\item{seed}{[integer, >0] initialization of the pseudo-random number generator.}

\item{keep.indiv}{[logical] should covariate specific occupancy probabilities be output?}

\item{trace}{[logical] should a progress bar be used to display the execution of the resampling procedure?}
}
\description{
Use multiple Cox models to fit an Illness Death Model from data in the wide format.
Output the fitted models and occupancy probabilities under different scenario (observed data, no transition to intermediate states, transition to a single intermediate state).
}
\examples{

library(survival)
library(mstate)
library(ggplot2)
library(riskRegression)

#### data (remove ties) ###
data(ebmt3)
set.seed(10)
noise <- sort(rnorm(NROW(ebmt3$prtime),sd = 0.00001))
ebmt3$prtime <- ebmt3$prtime/365.25 + noise
ebmt3$rfstime <- ebmt3$rfstime/365.25 + noise

#### PH without covariates ####
e.riskPH <- riskIDM(~1, data = ebmt3, PH = TRUE,
                    var.id = "id", n.boot = 100, seed = 10,
                    var.type = c("prstat", "rfsstat"),
                    var.time = c("prtime", "rfstime"))
summary(e.riskPH)
model.frame(e.riskPH)
model.tables(e.riskPH)
confint(e.riskPH, time = 1:10)
confint(e.riskPH, time = 2:11)
coef(e.riskPH, time = 1:10)
model.tables(e.riskPH, contrast = "all", time = c(1,2))
confint(e.riskPH, contrast = "all", time = c(1,2))
coef(e.riskPH, contrast = "all", time = c(1,2))
plot(e.riskPH, by = "scenario")
plot(e.riskPH, by = "scenario", scenario = "all")
plot(e.riskPH, by = "state")
plot(e.riskPH, by = "state", state = "all")
plot(e.riskPH, by = "contrast")

## parallel calculations
e.riskPH2 <- riskIDM(~1, data = ebmt3, PH = TRUE,
                    var.id = "id", n.boot = 100, cpus = 5, seed = 10,
                    var.type = c("prstat", "rfsstat"),
                    var.time = c("prtime", "rfstime"))

#### PH with covariates ####
dt.ebmt3 <- as.data.table(ebmt3)
eAdj.riskPH <- riskIDM(~age, data = dt.ebmt3, PH = TRUE,
                       var.id = "id", 
                       var.type = c("prstat", "rfsstat"),
                       var.time = c("prtime", "rfstime"))

coef(eAdj.riskPH, time = 1:10)
plot(eAdj.riskPH, by = "scenario")
plot(eAdj.riskPH, by = "state")
model.frame(eAdj.riskPH)

#### NPH without covariates ####
e.riskNPH <- riskIDM(~1, data = ebmt3, PH = FALSE,
                    var.id = "id", 
                    var.type = c("prstat", "rfsstat"),
                    var.time = c("prtime", "rfstime"))

model.frame(e.riskNPH)
plot(e.riskNPH)
plot(e.riskNPH, by = "state")
plot(e.riskNPH, by = "contrast")

## comparison with mstate
newdata.L <- data.frame(trans = c(1,2,3), strata = c(1,2,3))
tmat <- trans.illdeath(names = c("Tx", "PR", "RelDeath"))

msbmt <- msprep(time = c(NA, "prtime", "rfstime"),
                status = c(NA, "prstat", "rfsstat"),
                data = ebmt3, trans = tmat)

ebmt.coxNPH <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data = msbmt)
ebmt.coxPH <- coxph(Surv(Tstart, Tstop, status) ~ from + strata(to), data = msbmt)

ebmt.msfitNPH <- msfit(ebmt.coxNPH, newdata = newdata.L, trans = tmat)
ebmt.probNPH <- probtrans(ebmt.msfitNPH, predt = 0)
plot(ebmt.probNPH, use.ggplot = TRUE)

e.survNPH.obs <- confint(e.riskNPH, state = "survival")[scenario == "observed",]
plot(ebmt.probNPH[[1]]$time, ebmt.probNPH[[1]]$pstate1, type = "l")
points(e.survNPH.obs$time, e.survNPH.obs$estimate, type = "l", col = "red")

e.riskNPH.obs <- confint(e.riskNPH, state = "risk.2")[scenario == "observed",]
plot(ebmt.probNPH[[1]]$time, ebmt.probNPH[[1]]$pstate3, type = "l")
points(e.riskNPH.obs$time, e.riskNPH.obs$estimate, type = "l", col = "red")

#### NPH with covariates ####
eAdj.riskNPH <- riskIDM(~age, data = ebmt3, PH = FALSE,
                       var.id = "id", 
                       var.type = c("prstat", "rfsstat"),
                       var.time = c("prtime", "rfstime"))

summary(eAdj.riskNPH)
coef(eAdj.riskNPH)
plot(eAdj.riskNPH, by = "scenario")
plot(eAdj.riskNPH, by = "state", state = "all")

#### multiple exposures ####
set.seed(10)
n <- 1000 ## sample size (half)
tau <- 12 ## max follow-up time
Tevent <- rexp(2*n, rate = 1/10)
Tswitch <- c(runif(n, min = 1, max = 6), rep(Inf, n))
Cswitch <- sample.int(2, size = 2*n, replace = 2)
index.OC <- which(Tevent[1:n]>Tswitch[1:n])
Tevent[index.OC] <- Tswitch[index.OC] + rexp(length(index.OC), rate = c(1/5,1/2.5)[Cswitch[index.OC]])

df.W <- data.frame(id = 1:(2*n),
                   gender = 0:1,
                   time.event = pmin(Tevent,tau),
                   time.switch = pmin(Tevent,Tswitch,tau),
                   switch = ifelse(Tswitch<pmin(Tevent,tau), Cswitch, 0), 
                   event = as.numeric(Tevent <= tau))
df.W$switch <- factor(df.W$switch, levels = 0:2, c(0,"OC","IUD"))
# df.W$switch <- factor(df.W$switch, levels = 0:2, c("H","OC","IUD"))
# df.W$event <- factor(df.W$event, levels = 0:1, c("H","MDD"))
eME.riskPH <- riskIDM(~1, data = df.W, PH = FALSE,
                    var.id = "id",
                    var.type = c("switch","event"),
                    var.time = c("time.switch","time.event"))

summary(eME.riskPH)
model.frame(eME.riskPH)
confint(eME.riskPH, contrast = "all")
plot(eME.riskPH, by = "scenario")
plot(eME.riskPH, by = "scenario", scenario = "all")
plot(eME.riskPH, by = "state")
plot(eME.riskPH, by = "state", state = "all")
plot(eME.riskPH, by = "contrast", scenario = "all")
plot(eME.riskPH, by = "contrast", scenario = c("observed","no OC, IUD"))
plot(eME.riskPH, by = "contrast", scenario = c("observed","no OC, IUD","IUD instead of OC"))

#### exposure stratified on time of start ####

ebmt3$prstat3 <- factor(ebmt3$prstat * as.numeric(cut(ebmt3$prtime,c(0,0.06,0.1,10))))
table(ebmt3$prstat3, useNA = "always")

e.riskPH3 <- riskIDM(~1, data = ebmt3, PH = TRUE,
                    var.id = "id", n.boot = 100, seed = 10,
                    var.type = c("prstat3", "rfsstat"),
                    var.time = c("prtime", "rfstime"))

model.frame(e.riskPH3)
autoplot(e.riskPH3, by = "scenario") + coord_cartesian(xlim = c(0,0.2))
autoplot(e.riskPH3, stackplot = FALSE, by = "scenario") + coord_cartesian(xlim = c(0,0.2))


ebmt3.split <- rbind(data.frame(id = ebmt3$id, exposure = 0, start = 0,
                                 stop = pmin(ebmt3$rfstime, ebmt3$prtime),
                                 event = ebmt3$rfsstat*(ebmt3$prstat==0)),
                     data.frame(id = ebmt3[ebmt3$prstat==1,"id"], exposure = 1,
                                 start = ebmt3[ebmt3$prstat==1,"prtime"],
                                 stop = ebmt3[ebmt3$prstat==1,"rfstime"],
                                 event = ebmt3[ebmt3$prstat==1,"rfsstat"])
)

model.frame(e.riskPH)
coxph(Surv(start,stop,event)~exposure, data = ebmt3.split)

model.frame(e.riskPH3)
coxph(Surv(start,stop,event)~exposure, data = ebmt3.split[ebmt3.split$start<=0.06,])
}
