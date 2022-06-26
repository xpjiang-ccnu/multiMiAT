\name{multiMiAT}
\alias{multiMiAT}
\title{multiMiAT test}

\usage{
multiMiAT(y, otu.tab, tree, covs = NULL, model = c("bcl", "clm"), distance = c("BC", "U.UniFrac", "G.UniFrac.00", "G.UniFrac.25", "G.UniFrac.50", "G.UniFrac.75", "W.UniFrac"), kernel = c("LK", "GK", "LaK"), W, n.perm = 999, seed = 1)
}

\description{
This function can be used for testing the association between microbiome and multicategory phenotypes. Specifically, multiMiAT can be applied to not only ordinal multicategory phenotypes (i.e., disease severity), but also nominal multicategory phenotypes (i.e., tumor subtype or dietary pattern). Considering that the binary outcome is a special case of the multicategory outcomes, our method can also be used for binary outcomes.
}

\arguments{
  \item{y}{
Multicategory outcomes (i.e., host phenotype of interest). Multicategory outcomes mainly contain ordinal multicategory outcomes (e.g., disease severity) nomial multicategory outcomes (e.g., tumor subtype).
}
  \item{otu.tab}{
OTU count table. Rows are samples and columns are OTUs.
}
  \item{tree}{
Phylogenetic tree. If NULL, the phylogeny-based distance cannot be calculated.
}
  \item{covs}{
Confounding factors (e.g., age, gender). Default is covs = NULL.
}
  \item{model}{
Multinomial logit models. Default is model = c("bcl", "clm"). "bcl" for baseline category logit model, "clm" for cumulative link model. In addition, "crm" for continuation ratio model can also be used.
}
  \item{distance}{
Distance. Default is distance = c("BC", "U.UniFrac", "G.UniFrac.00", "G.UniFrac.25", "G.UniFrac.50", "G.UniFrac.75", "W.UniFrac"). "BC" for Bray-Curtis distance, "U.UniFrac" for Unifrac distance, "G.UniFrac.00" for Generalized UniFrac distance with alpha = 0, "G.UniFrac.25" for Generalized UniFrac distance with alpha = 0.25, "G.UniFrac.50" for Generalized UniFrac distance with alpha = 0.5, "G.UniFrac.75" for Generalized UniFrac distance with alpha = 0.75, "W.UniFrac" for Weighted UniFrac distance.
}
  \item{kernel}{
Kernel function. Default is kernel = c("LK", "GK", "LaK"). "LK" for linear kernel, "GK" for Gaussian kernel, "LaK" for Laplacian kernel.
}
  \item{W}{
Weight vector of OTUs.
}
  \item{n.perm}{
A number of permutations. Default is n.perm = 999. 
}
  \item{seed}{
Random number generator. Default is seed = 1.
}
}

\value{
multiMiRKAT.bcl.pvs: The p-values for the individual tests based on baseline category logit model.

multiMiRKAT.clm.pvs: The p-values for the individual tests based on cumulative link model.

Score.test: The p-value for the bcl-based score test.

MiRKAT.MC: The p-values for the MiRKAT-MC tests (i.e., MiRKAT-MCN and MiRKAT-MCO).

multiMiRKAT: The p-values for the multiMiRKAT tests (i.e., multiMiRKAT-N and multiMiRKAT-O).

multiMiAT: The p-value for the optimal test (i.e., multiMiAT).
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
W <- matrix(1, nrow = 1, ncol = ncol(otu.tab))

# Fit multiMiAT
\dontrun{
out <- multiMiAT(y = y, otu.tab = otu.tab, tree = tree, covs = covs, model = c("bcl", "clm"), distance = c("BC", "U.UniFrac", "G.UniFrac.00", "G.UniFrac.25", "G.UniFrac.50", "G.UniFrac.75", "W.UniFrac"), kernel = c("LK", "GK", "LaK"), W = W)}
}
