\name{multiMiRKAT_bcl}
\alias{multiMiRKAT_bcl}
\title{multiMiRKAT based on baseline category logit model}

\usage{
multiMiRKAT_bcl(y, Ks, covs = NULL, n.perm = 999, seed = 1)
}

\description{
This function can be used to detect the association between the microbiome and multinomial outcomes, which is microbiome regression-based kernel association tests based on baseline category logit model.
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
  \item{n.perm}{
A number of permutations. Default is n.perm = 999. 
}
  \item{seed}{
Random number generator. Default is seed = 1.
}
}

\value{
$T.statistic: The test statistic of the individual tests based on baseline category logit model.

$multiMiRKAT.bcl.pvs: The p-values for the individual tests based on baseline category logit model.

$multiMiRKAT.bcl.Omnibus.p: The p-values for the omnibus tests based on baseline category logit model.
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

# Fit multiMiRKAT_bcl
\dontrun{
Ks <- Whole.Kernels(otu.tab, tree) 
out <- multiMiRKAT_bcl(y = y, Ks = Ks, covs = covs, distance = c("BC", "U.UniFrac", "G.UniFrac.00", "G.UniFrac.25", "G.UniFrac.50", "G.UniFrac.75", "W.UniFrac"), kernel = c("LK", "GK", "LaK"))}
}
