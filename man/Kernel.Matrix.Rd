\name{Kernel.Matrix}
\alias{Kernel.Matrix}
\title{Kernel matrix obtained through the kernel function}

\usage{
Kernel.Matrix(D, kernel, sigma.GK, sigma.EK, sigma.LaK)
}

\description{
This function establishes the kernel matrix through the kernel function (i.e., linear kernel, Gaussian kernel, exponential kernel and Laplacian kernel function).
}

\arguments{
  \item{D}{
Distance matrix.
}
  \item{kernel}{
Kernel function. Default is kernel = c("LK", "GK", "EK", "LaK"). "LK" for linear kernel, "GK" for Gaussian kernel, "EK" for exponential kernel, "LaK" for Laplacian kernel.
}
  \item{sigma.GK}{
Parameter setting of Gaussian kernel function.
}
  \item{sigma.EK}{
Parameter setting of exponential kernel function.
}
  \item{sigma.LaK}{
Parameter setting of Laplacian kernel function.
}
}

\value{
$LK: Kernel matrix using linear kernel.

$GK.sigma1: Kernel matrix using Gaussian kernel with sigma = 1.

$EK.sigma1: Kernel matrix using exponential kernel with sigma = 1.

$LaK.sigma1: Kernel matrix using Laplacian kernel with sigma = 1.
}

\examples{
require(GUniFrac)
data(throat.otu.tab)
data(throat.tree)
unifs <- GUniFrac::GUniFrac(throat.otu.tab, throat.tree, alpha = c(0, 0.25, 0.5, 0.75, 1))$unifracs
u.unif <- unifs[, , "d_UW"]
\dontrun{
u.unif.k <- Kernel.Matrix(u.unif, kernel = c("LK", "GK", "EK", "LaK"), sigma.GK = 1, sigma.EK = 1, sigma.LaK = 1)}
}
