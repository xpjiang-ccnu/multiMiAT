\name{score.test_bcl}
\alias{score.test_bcl}
\title{bcl-based score test}

\usage{
score.test_bcl(y, otu.tab, covs, W)
}

\description{
This function is score test based on baseline category logit model, which is a special case of an individual test in subtype analysis with somatic mutations (SASOM). It can be also used for the association between the microbiome and multicategory outcomes.
}

\arguments{
  \item{y}{
Multicategory outcomes (i.e., host phenotype of interest). Multicategory outcomes mainly contain ordinal multicategory outcomes (e.g., disease severity) nomial multicategory outcomes (e.g., tumor subtype). The data format can be factor, numeric, integer, et al.
}
  \item{otu.tab}{
OTU count table. Rows are samples and columns are OTUs.
}
  \item{covs}{
Confounding factors (e.g., age, gender). Default is covs = NULL. The data format must be data frame.
}
  \item{W}{
Weight vector of OTUs.
}
}

\value{
The p-value for the score test based on baseline category logit model.
}

\examples{
# Import requisite R packages
require(phyloseq)
require(VGAM)
require(MASS)
require(permute)
require(ape)
require(Matrix)

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

# Fit bcl-based score test
\dontrun{
out <- score.test_bcl(y, otu.tab, covs, W)}
}
