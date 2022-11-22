\name{MiRKAT_MCO}
\alias{MiRKAT_MCO}
\title{MiRKAT-MCO test}

\usage{
out.MiRKAT_MCO <- MiRKAT_MCO(y, Ks, covs)
}

\description{
MiRKAT-MCO test is MiRKAT-MC test based on cumulative link model, which mainly refers to the MiRKATMC function that sets "data.type" to "ordinal" in R package MiRKATMC. In addition, we fix a bug in the MiRKAT function, i.e., the p-value combination method HMP reports an error when zero value appears in these combined p-values.
}

\arguments{
  \item{y}{
Multicategory outcomes (i.e., host phenotype of interest). Multicategory outcomes mainly contain ordinal multicategory outcomes (e.g., disease severity) nomial multicategory outcomes (e.g., tumor subtype). The data format can be factor, numeric, integer, et al.
}
  \item{Ks}{
Kernel matrix based on diverse distances and kernel functions.
}
  \item{covs}{
Confounding factors (e.g., age, gender). Default is covs = NULL. The data format must be data frame.
}
}

\value{
MiRKAT_MCO.pvs: The p-values for the individual tests of MiRKAT-MCO test.

MiRKAT_MCO: The p-value for the omnibus test of MiRKAT-MCO test.
}

\examples{
# Import requisite R packages
require(ecodist)
require(GUniFrac)
require(phyloseq)
require(VGAM)
require(MASS)
require(permute)
require(ape)
require(Matrix)
require(multiMiAT)

# Import example microbiome data
data(biom1)
otu.tab <- otu_table(biom1)
tree <- phy_tree(biom1)
sample.data <- sample_data(biom1)
y <- sample.data$DiseaseState
y[which(y == "Health")] = 1
y[which(y == "Adenoma")] = 2
y[which(y == "Cancer")] = 3
y <- as.factor(y)
x1 <- sample.data$Age
x2 <- sample.data$Race
x2[which(x2 == "White")] <- 1
x2[which(x2 == "Other")] <- 0
x2 <- as.factor(x2)
covs <- as.data.frame(cbind(x1, x2))
colnames(covs) <- c("Age", "Race")
covs[,1] <- as.numeric(covs[,1])
covs[,2] <- as.factor(covs[,2])

# Fit MiRKAT-MCO
\dontrun{
Ks <- Whole.Kernels(otu.tab, tree, distance = c("BC", "U.UniFrac", "G.UniFrac.00", "G.UniFrac.25", "G.UniFrac.50", "G.UniFrac.75", "W.UniFrac"), kernel = c("LK", "GK", "LaK"))
out.MiRKAT_MCO <- MiRKAT_MCO(y = y, Ks = Ks, covs = covs)}
}
