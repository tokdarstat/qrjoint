\name{chull.center}
\Rdversion{0.1-1}
\alias{chull.center}
\title{Fast Interior Point Center of Multivariate Data}
\description{Calculates an interior point by averaging a small number of near-extreme points of the cloud.} 
\usage{ 
chull.center(x, maxEPts = ncol(x) + 1, plot = FALSE)
} 
\arguments{
\item{x}{a matrix giving the data cloud.}
\item{maxEPts}{integer giving the maximum number of (near)-extreme points to be used in averaging. Default is \code{ncol(x)+1}.}
\item{plot}{logical indicating whether a pairwise scatter plot should be made}
}
\value{
Returns an interior point of the data cloud. The positions of the near extreme points are returned as the attribute \code{"EPts"}.
}

\details{Near extreme points are found in a space-filling manner by adding points with minimum residual conditional variance given points already included under a smooth GP specification. See Yang and Tokdar (2015), Section B.1. for more details.}


\examples{ 
p <- 9
n <- 200
u <- runif(n)
require(splines)
x <- bs(u, df = p)
chull.center(x, plot = TRUE)
}
\keyword{programming}