\name{Whole.Kernels}
\alias{Whole.Kernels}
\title{Construction of the kernel matrix through diverse distances}

\usage{
Whole.Kernels(otu.tab, tree, distance = c("BC", "U.UniFrac", "G.UniFrac.00", "G.UniFrac.25", "G.UniFrac.50", "G.UniFrac.75", "W.UniFrac"), kernel = c("LK", "GK", "LaK"), sigma.GK = 1, sigma.EK = 1, sigma.LaK = 1)
}

\description{
The construction of the kernel matrix through diverse distances (i.e., Jaccard dissimilarity, Bray-Curtis distance, Unifrac distance, Generalized UniFrac distance with alpha = 0, Generalized UniFrac distance with alpha = 0.25, Generalized UniFrac distance with alpha = 0.5, Generalized UniFrac distance with alpha = 0.75 and Weighted UniFrac distance) and kernel functions (i.e., linear kernel, Gaussian kernel, exponential kernel and Laplacian kernel function).
}

\arguments{
  \item{otu.tab}{
OTU count table. Rows are samples and columns are OTUs.
}
  \item{tree}{
Phylogenetic tree. If NULL, the phylogeny-based distance cannot be calculated.
}
  \item{distance}{
Distance. Default is distance = c("BC", "U.UniFrac", "G.UniFrac.00", "G.UniFrac.25", "G.UniFrac.50", "G.UniFrac.75", "W.UniFrac"), kernel = c("LK", "GK", "EK", "LaK"). "Jaccard" for Jaccard dissimilarity, "BC" for Bray-Curtis distance, "U.UniFrac" for Unifrac distance, "G.UniFrac.00" for Generalized UniFrac distance with alpha = 0, "G.UniFrac.25" for Generalized UniFrac distance with alpha = 0.25, "G.UniFrac.50" for Generalized UniFrac distance with alpha = 0.5, "G.UniFrac.75" for Generalized UniFrac distance with alpha = 0.75, "W.UniFrac" for Weighted UniFrac distance. In addition, "Jaccard" is also alternative distance.
}
  \item{kernel}{
Kernel function. Default is kernel = c("LK", "GK", "LaK"). "LK" for linear kernel, "GK" for Gaussian kernel, "EK" for exponential kernel, "LaK" for Laplacian kernel. In addition, "EK" is also alternative kernel function.
}
  \item{sigma.GK}{
Parameter setting of Gaussian kernel function. Default is sigma.GK = 1.
}
  \item{sigma.EK}{
Parameter setting of exponential kernel function. Default is sigma.EK = 1.
}
  \item{sigma.LaK}{
Parameter setting of Laplacian kernel function. Default is sigma.LaK = 1.
}
}

\value{
$Jaccard: Kernel matrix based on Jaccard dissimilarity and diverse kernel functions (i.e., linear kernel, Gaussian kernel, exponential kernel and Laplacian kernel function).

$BC: Kernel matrix based on Bray-Curtis distance and diverse kernel functions (i.e., linear kernel, Gaussian kernel, exponential kernel and Laplacian kernel function).

$U.UniFrac: Kernel matrix based on Unifrac distance and diverse kernel functions (i.e., linear kernel, Gaussian kernel, exponential kernel and Laplacian kernel function).

$G.UniFrac.00: Kernel matrix based on Generalized UniFrac distance with alpha = 0 and diverse kernel functions (i.e., linear kernel, Gaussian kernel, exponential kernel and Laplacian kernel function).

$G.UniFrac.25: Kernel matrix based on Generalized UniFrac distance with alpha = 0.25 and diverse kernel functions (i.e., linear kernel, Gaussian kernel, exponential kernel and Laplacian kernel function).

$G.UniFrac.50: Kernel matrix based on Generalized UniFrac distance with alpha = 0.5 and diverse kernel functions (i.e., linear kernel, Gaussian kernel, exponential kernel and Laplacian kernel function).

$G.UniFrac.75: Kernel matrix based on Generalized UniFrac distance with alpha = 0.75 and diverse kernel functions (i.e., linear kernel, Gaussian kernel, exponential kernel and Laplacian kernel function).

$W.UniFrac: Kernel matrix based on Weighted UniFrac distance and diverse kernel functions (i.e., linear kernel, Gaussian kernel, exponential kernel and Laplacian kernel function).
}

\examples{
require(GUniFrac)
require(ecodist)
data(throat.otu.tab)
data(throat.tree)
\dontrun{
Ks <- Whole.Kernels(throat.otu.tab, throat.tree)}
}
