\name{ordinal.simulated}
\alias{ordinal.simulated}
\title{Simulated data including OTU count table, ordinal multinomial host phenotypes and confounding factors}

\usage{
ordinal.simulated(betas, proportion, xformula, xdata, link)
}

\description{
This function generates the simulated data including OTU count table, ordinal multinomial host phenotypes (e.g., disease severity) and confounding factors (e.g., age, gender).
}

\arguments{
  \item{betas}{
Marginal regression parameter vector.
}
  \item{proportion}{
Proportion of each category. Balanced and unbalanced experiment data can be designed by setting this parameter.
}
  \item{xformula}{
Formula expression without including outcomes.
}
  \item{xdata}{
Data containing the variables provided in xformula.
}
  \item{link}{
Link function of multinomial logit models. Default is link = "logit".
}
}

\value{
@otu\_table: OTU counts table simulated based on real data.

@tax\_table: Taxa information.

@sam\_data: Simulated clinical information, including multinomial host phenotypes of interest and confounding factors.

@phy\_tree: Phylogenetic tree.

@refseq: A biological sequence set object of a class.
}

\examples{
# Parameter settings:
sample_size <- 100
OTU_select <- 30
betas0 <- c(0.5, 0.5)
beta <- 1

# Generation of microbiome data:
require(phyloseq)
data("phy", package = "MiHC")
phy.otu.tab <- otu_table(phy)
tree <- phy_tree(phy)
# library(dirmult)
# phy.parameters <- dirmult(phy.otu.tab)
data("phy.parameters")
otu.tab <- SimulateOTU(phy.otu.tab, nSam = 100, phy.parameters, mu = 1000, size = 25)

# Random selection of OTUs:
otu.tab_selection <- otu.tab[, order(apply(otu.tab, 2, sum), decreasing = TRUE)[1:OTU_select]]
colnames(otu.tab_selection) <- paste("OTU", colnames(otu.tab_selection), sep = "")
otu.tab_selection <- scale(otu.tab_selection)

# Generation of covariates and regression parameters:
betas1 <- rep(beta, ncol(otu.tab_selection))
x1 <- rbinom(sample_size, 1, 0.5)
x2 <- rnorm(sample_size, 0, 1)
x1 <- scale(x1)
x2 <- scale(x2)
xdata <- data.frame(x1, x2, otu.tab_selection)
beta_coefficients <- c(betas0, betas1)
xnam <- colnames(xdata)
xformula <- as.formula(paste("~ ", paste(xnam, collapse= "+")))
covs <- data.frame(x1, x2)
covs[,1] <- as.factor(covs$x1) 
covs[,2] <- as.numeric(covs$x2)

# Generation of ordinal multinomial outcomes:
y <- ordinal.simulated(betas = beta_coefficients, proportion = c(1, 1, 1, 1), xformula = xformula, xdata = xdata, link = "logit") 

# Data synthesis of simulated data with ordinal multinomial outcome:
otu.tab_information <- otu_table(otu.tab, taxa_are_rows = FALSE)
sample_information <- sample_data(as.data.frame(cbind(y, x1, x2)))
colnames(sample_information) <- c("y", "x1", "x2")
rownames(sample_information) <- rownames(otu.tab_information)
biom <- phyloseq(otu.tab_information, tree, sample_information)
}
