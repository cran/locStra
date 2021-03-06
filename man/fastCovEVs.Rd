% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/locStra.r
\name{fastCovEVs}
\alias{fastCovEVs}
\title{Computation of the k leading eigenvectors of the covariance matrix for a (sparse) input matrix.}
\usage{
fastCovEVs(m, k, useCpp = TRUE, sparse = TRUE, q = 2)
}
\arguments{
\item{m}{A (sparse) matrix for which the eigenvectors of its covariance matrix are sought. The input matrix is assumed to be oriented to contain the data for one individual per column.}

\item{k}{The number of leading eigenvectors.}

\item{useCpp}{Flag to switch between R or C++ implementations. Default is \code{useCpp=TRUE}.}

\item{sparse}{Flag to switch between purpose-built dense or sparse implementations. Default is \code{sparse=TRUE}.}

\item{q}{The number of power iteration steps (default is \code{q=2}).}
}
\value{
The k leading eigenvectors of the covariance matrix of \code{m} as a column matrix.
}
\description{
Computation of the k leading eigenvectors of the covariance matrix for a (sparse) input matrix.
}
\examples{
require(locStra)
require(Matrix)
m <- matrix(sample(0:1,100,replace=TRUE),ncol=5)
sparseM <- Matrix(m,sparse=TRUE)
print(fastCovEVs(sparseM,k=2,useCpp=FALSE))

}
\references{
R Core Team (2014). R: A Language and Environment for Statistical Computing. R Foundation for Stat Comp, Vienna, Austria.

N. Halko, P.G. Martinsson, and J.A. Tropp (2011). Finding Structure with Randomness: Probabilistic Algorithms for Constructing Approximate Matrix Decompositions. SIAM Review: 53(2), pp. 217--288.
}
