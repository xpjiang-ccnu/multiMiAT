\name{sigma_value}
\alias{sigma_value}
\title{Median value of distance matrix square or distance matrix}

\usage{
sigma_value(D)
}

\description{
The median value of distance matrix square (D * D) or distance matrix (D).
}

\arguments{
  \item{D}{
Distance matrix.
}
}

\value{
$sigma.GK: Parameter sigma setting of Gaussian kernel function.

$sigma.EK: Parameter sigma setting of exponential kernel function.

$sigma.LaK: Parameter sigma setting of Laplacian kernel function.
}

\examples{
require(GUniFrac)
data(throat.otu.tab)
data(throat.tree)
unifs <- GUniFrac::GUniFrac(throat.otu.tab, throat.tree, alpha = c(0, 0.25, 0.5, 0.75, 1))$unifracs
u.unif <- unifs[, , "d_UW"]
sigma.u.unif.list <- sigma_value(u.unif)
}
