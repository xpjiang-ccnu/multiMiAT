\name{SimulateOTU}
\alias{SimulateOTU}
\title{OTUs count table simulated based on the Dirichlet-multinomial model}

\usage{
SimulateOTU(data, nSam, parameters, mu, size)
}

\description{
This function generates the OTU count table simulated based on the Dirichlet-multinomial model according to real data, which is the same as the SimulateOTU function of the GEEMiHC package. We first calculate the parameters of the real data, and then generate the simulated data according to the parameters.
}

\arguments{
  \item{data}{
Real data.
}
  \item{nSam}{
Sample size.
}
  \item{parameters}{
The estimated parameter based on a real microbiome data, including OTU proportions and overdispersion parameter.

}
  \item{mu}{
The mean of the negative binomial distribution.
}
  \item{size}{
The size of the negative binomial distribution.
}
}

\value{
OTU: OTU counts table simulated based on real data.
}

\examples{
require(phyloseq)
data("phy", package = "MiHC")
phy.otu.tab <- otu_table(phy)
# require(dirmult)
# phy.parameters <- dirmult(phy.otu.tab)
data("phy.parameters")
otu.tab <- SimulateOTU(phy.otu.tab, nSam = 100, phy.parameters, mu = 1000, size = 25)
}
