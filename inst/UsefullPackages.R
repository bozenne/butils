packages <- c(
  "car","caret", "cmprsk", "coda",
  "data.table","devtools","doMC",
  "e1071",
  "fields","foreign",
  "glmnet",
  "MASS","Matrix","mgcv",
  "geepack","ggmap","ggplot2","gof",
  "knitr",
  "lava","lava.tobit","lme4",
  "maps","maptools","markdown","mritc","multcomp","mvtnorm",
  "nlme","nnet"," 	numDeriv",
  "parallel","pracma",
  "optimx","oro.dicom", "oro.nifti",
  "party","pec","penalized","proc","prodlim","psych","publish",
  "quantreg",
  "RANN","rbenchmark","Rcpp","RcppArmadillo","RcppEigen","ReporteRs","rgl","riskRegression","rms,","ROCR","roxygen2","rpart",
  "shiny","snowfall","sp","spam","stringr",
  "testthat","tcltk","timereg",
  "ranger","randomForest","randomForestSRC",
  "survival",
  "xlsx","xtable"
)


http://askubuntu.com/questions/614530/how-to-install-latest-version-of-r-on-ubuntu-12-04-lts
codename=$(lsb_release -c -s)
echo "deb http://cran.fhcrc.org/bin/linux/ubuntu $codename/" | sudo tee -a /etc/apt/sources.list > /dev/null
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
sudo add-apt-repository ppa:marutter/rdev
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install r-base r-base-dev

