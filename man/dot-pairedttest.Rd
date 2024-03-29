% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pairedCompMean.R
\name{.pairedttest}
\alias{.pairedttest}
\title{Fast Paired t-test}
\usage{
.pairedttest(x, y)
}
\arguments{
\item{x}{vector of numeric values}

\item{y}{vector of numeric values (reference)}
}
\description{
Fast paired t-test
}
\examples{
X <- rnorm(1e2)
Y <- rnorm(1e2)
t.test(X,Y,paired = TRUE)
.pairedttest(X,Y)
microbenchmark::microbenchmark(t.test(X,Y,paired = TRUE), .pairedttest(X,Y))
}
\keyword{internal}
