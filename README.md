# multiMiAT

Type: Package

Title: multiMiAT: An optimal microbiome-based association test for multicategory phenotypes

Version: 1.0

Author: Han Sun

Maintainer: Han Sun sunh529@mails.ccnu.edu.cn; Xingpeng Jiang xpjiang@mail.ccnu.edu.cn

Imports: phyloseq, ecodist, GUniFrac, permute, ape, VGAM, MASS, Matrix, dirmult, harmonicmeanp, devtools, MiHC

Description: multiMiAT can be used to test the association between microbiome and multicategory phenotypes.

License: GPL-2

Encoding: UTF-8

LazyData: true

URL: https://github.com/xpjiang-ccnu/multiMiAT



## Introduction
This R package, **multiMiAT**, can be used for testing the association between microbiome and multicategory phenotypes.  Specifically, **multiMiAT** can be applied to not only ordinal multicategory phenotypes (i.e., disease severity), but also nominal multicategory phenotypes (i.e., tumor subtype or dietary pattern). Considering that the binary outcome is a special case of the multicategory outcomes, our method can also be used for binary outcomes.



## Installation 
phyloseq:
```
BiocManager::install("phyloseq")
```

ecodist:
```
install.packages("ecodist")
```

GUniFrac:
```
install.packages("GUniFrac")
```

permute:
```
install.packages("permute")
```

ape:
```
install.packages("ape")
```

VGAM:
```
install.packages("VGAM")
```

MASS:
```
install.packages("MASS")
```

Matrix:
```
install.packages("Matrix")
```

dirmult:
```
install.packages("dirmult")
```

harmonicmeanp:
```
install.packages("harmonicmeanp")
```

devtools:
```
install.packages("devtools")
```

MiHC:
```
devtools::install_github('hk1785/MiHC', force = T)
```


You may install `multiMiAT` from GitHub using the following code: 

```
devtools::install_github("xpjiang-ccnu/multiMiAT", force=T)
```

---------------------------------------------------------------------------------------------------------------------------------------



## Usage
```
multiMiAT(y, otu.tab, tree, covs = NULL, method = c("multiMiRKAT", "score", "MiRKAT-MC"), distance = c("BC", "U.UniFrac", "G.UniFrac.00", "G.UniFrac.25", "G.UniFrac.50", "G.UniFrac.75", "W.UniFrac"), kernel = c("LK", "GK", "LaK"), W, n.perm = 999, seed = 1) 
```


## Arguments
* _y_ - Multicategory outcomes (i.e., host phenotype of interest). Multicategory outcomes mainly contain ordinal multicategory outcomes (e.g., disease severity) nomial multicategory outcomes (e.g., tumor subtype). The data format can be factor, numeric, integer, et al.

* _otu.tab_ - OTU count table. Rows are samples and columns are OTUs.

* _tree_ - Phylogenetic tree. If NULL, the phylogeny-based distance cannot be calculated.

* _covs_ - Confounding factors (e.g., age, gender). Default is covs = NULL. The data format must be data frame.

* _method_ - The combined three methods (i.e., multiMiRKAT, score test and MiRKAT-MC). Any combination of the three tests can be performed. Default is method = c("multiMiRKAT", "score", "MiRKAT-MC") for multiMiAT test.

* _distance_ - Distance. Default is distance = c("BC", "U.UniFrac", "G.UniFrac.00", "G.UniFrac.25", "G.UniFrac.50", "G.UniFrac.75", "W.UniFrac"). "BC" for Bray-Curtis distance, "U.UniFrac" for Unifrac distance, "G.UniFrac.00" for Generalized UniFrac distance with alpha = 0, "G.UniFrac.25" for Generalized UniFrac distance with alpha = 0.25, "G.UniFrac.50" for Generalized UniFrac distance with alpha = 0.5, "G.UniFrac.75" for Generalized UniFrac distance with alpha = 0.75, "W.UniFrac" for Weighted UniFrac distance.

* _kernel_ - Kernel function. Default is kernel = c("LK", "GK", "LaK"). "LK" for linear kernel, "GK" for Gaussian kernel, "LaK" for Laplacian kernel.

* _W_ - Weight vector of OTUs.

* _CLR_ - Centered log-ratio (CLR) transformation. Default is CLR = FALSE for no CLR transformation.

* _n.perm_ - A number of permutations. Default is n.perm = 999. 

* _seed_ - Random number generator. Default is seed = 1.


## Values
_$multiMiRKAT.pvs_ - The _p_-values for the individual tests of multiMiRKAT.

_$Combined.methods.pvs_ - The _p_-values for the combined tests (i.e., multiMiRKAT, score test and MiRKAT-MC).

_$Optimal.test.p_ - The _p_-value for the optimal test (i.e., multiMiAT).


## Example
**Import requisite R packages:**

```
require(ecodist)
require(GUniFrac)
require(phyloseq)
require(VGAM)
require(MASS)
require(permute)
require(ape)
require(Matrix)
require(multiMiAT)
```

**Import example microbiome data:**

```
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
```

**Fit multiMiAT:**

```
out <- multiMiAT(y = y, otu.tab = otu.tab, tree = tree, covs = covs, method = c("multiMiRKAT", "score", "MiRKAT-MC"), distance = c("BC", "U.UniFrac", "G.UniFrac.00", "G.UniFrac.25", "G.UniFrac.50", "G.UniFrac.75", "W.UniFrac"), kernel = c("LK", "GK", "LaK"), W = W)
out
```



## References
* Sun H, et al. multiMiAT: An optimal microbiome-based association test for multicategory phenotypes. (under review)

* Agarwal D, Zhang NR. Semblance: An empirical similarity kernel on probability spaces. _Science Advances_ 2019;**5**(12):eaau9630.

* Jiang Z, et al. MiRKAT-MC: A distance-based microbiome kernel association test with multi-categorical outcomes. _Frontiers in Genetics_ 2022;**13**:841764.

* Koh H, et al. A powerful microbiome-based association test and a microbial taxa discovery framework for comprehensive association mapping. _Microbiome_ 2017;**5**(1):45.

* Liu M, et al. A method for subtype analysis with somatic mutations. _Bioinformatics_ 2021;**37**(1):50-56.

* Wilson DJ. The harmonic mean p-value for combining dependent tests. _Proceedings of the National Academy of Sciences_ 2019;**116**(4):1195–1200.

* Yee TW, et al. The VGAM package for categorical data analysis. _Journal of Statistical Software_ 2010;**32**(10).

* Zhao N, et al. Testing in microbiome-profiling studies with MiRKAT, the microbiome regression-based kernel association test. _American Journal of Human Genetics_ 2015;**96**(5):797-807.



---------------------------------------------------------------------------------------------------------------------------------------



# Other functions

## Function **SimulateOTU**

### Description
The function, **SimulateOTU**, generates the OTU count table simulated based on the Dirichlet-multinomial model according to real data, which is the same as the **SimulateOTU** function of the **GEEMiHC** package. We first calculate the parameters of the real data, and then generate the simulated data according to the parameters.


### Usage
```
SimulateOTU(data, nSam, parameters, mu, size)
```


### Arguments
* _data_ - real data.

* _nSam_ - Sample size.

* _parameters_ - The estimated parameter based on a real microbiome data, including OTU proportions and overdispersion parameter.

* _mu_ - The mean of the negative binomial distribution.

* _size_ - The size of the negative binomial distribution.


### Values
_$OTU_ - OTU counts table simulated based on real data.


### Example
**Parameter estimation:**

```
require(phyloseq)
require(multiMiAT)
data("phy", package = "MiHC")
phy.otu.tab <- otu_table(phy)
data("phy.parameters")
```

**Generation of microbiome data:**

```
otu.tab <- SimulateOTU(phy.otu.tab, nSam = 100, phy.parameters, mu = 1000, size = 25)
```


### References
* Chen J and Li H. Variable selection for sparse Dirichlet-multinomial regression with an application to microbiome data analysis. _Annals of Applied Statistics_ 2013;**7**(1).

* Koh H and Zhao N. A powerful microbial group association test based on the higher criticism analysis for sparse microbial association signals. _Microbiome_ 2020;**8**(1):63.

* Sun H, et al. Detecting sparse microbial association signals adaptively from longitudinal microbiome data based on generalized estimating equations. _Briefings in Bioinformatics_ 2022:bbac149.

* Wu C, et al. An adaptive association test for microbiome data. _Genome Med_ 2016;**8**(1):56.



## Function **ordinal.simulated**

### Description
The function, **ordinal.simulated**, generates the simulated data including OTU count table, ordinal multinomial host phenotype (e.g., disease severity) and confounding factors (e.g., age, gender), which mainly refers to the **rmult.clm** function in R package, **SimCorMultRes**.


### Usage
```
ordinal.simulated(betas, proportion, xformula, xdata, link)
```


### Arguments
* _betas_ - Marginal regression parameter vector.

* _proportion_ - Proportion of each category. Balanced and unbalanced experiment data can be designed by setting this parameter.

* _xformula_ - Formula expression without including outcomes.

* _xdata_ - Data containing the variables provided in xformula.

* _link_ - Link function of multinomial logit models. Default is link = "logit".


### Values
_@otu\_table_ - OTU counts table simulated based on real data.

_@tax\_table_ - Taxonomy information.

_@sam\_data_ - Simulated clinical information, including multinomial host phenotypes of interest and confounding factors.

_@phy\_tree_ - Phylogenetic tree.

_@refseq_ - A biological sequence set object of a class.


### Example
**Parameter settings:**

```
sample_size <- 100
OTU_select <- 30
betas0 <- c(0.5, 0.5)
beta <- 1
```

**Generation of microbiome data:**

```
require(phyloseq)
require(multiMiAT)
data("phy", package = "MiHC")
phy.otu.tab <- otu_table(phy)
tree <- phy_tree(phy)
data("phy.parameters")
otu.tab <- SimulateOTU(phy.otu.tab, nSam = 100, phy.parameters, mu = 1000, size = 25)
```

**Random selection of OTUs:**

```
otu.tab_selection <- otu.tab[, order(apply(otu.tab, 2, sum), decreasing = T)[1:OTU_select]]
colnames(otu.tab_selection) <- paste("OTU", colnames(otu.tab_selection), sep = "")
otu.tab_selection <- scale(otu.tab_selection)
```

**Generation of covariates and regression parameters:**

```
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
```

**Generation of ordinal multicategory outcomes:**

```
y <- ordinal.simulated(betas = beta_coefficients, proportion = c(1, 1, 1, 1), xformula = xformula, xdata = xdata, link = "logit") 
```

**Data synthesis of simulated data with ordinal multicategory outcomes:**

```
otu.tab_information <- otu_table(otu.tab, taxa_are_rows = F)
sample_information <- sample_data(as.data.frame(cbind(y, x1, x2)))
colnames(sample_information) <- c("y", "x1", "x2")
rownames(sample_information) <- rownames(otu.tab_information)
biom <- phyloseq(otu.tab_information, tree, sample_information)
```


### References
* Bi Wenjian, et al. Efficient mixed model approach for large-scale genome-wide association studies of ordinal categorical phenotypes. _American Journal of Human Genetics_ 2021;**108**(5):825–839. 

* McMurdie PJ and Holmes S. phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. PLoS ONE. 2013;**8**(4):e61217.

* Touloumis A. Simulating correlated binary and multinomial responses under marginal model specification: The SimCorMultRes package. _The R Journal_ 2016,**8**(2):79.



## Function **nominal.simulated**

### Description
The function, **nominal.simulated**, generates the simulated data including OTU count table, nominal multinomial host phenotype (e.g., tumor subtype) and confounding factors (e.g., age, gender), which mainly refers to the **rmult.bcl** function in R package, **SimCorMultRes**.


### Usage
```
nominal.simulated(betas, categories, xformula, xdata)
```


### Arguments
* _betas_ - Marginal regression parameter vector.

* _categories_ - Number of categories.

* _xformula_ - Formula expression without including outcomes.

* _xdata_ - Data containing the variables provided in xformula.


### Values
_@otu\_table_ - OTU counts table simulated based on real data.

_@tax\_table_ - Taxa information.

_@sam\_data_ - Simulated clinical information, including multinomial host phenotypes of interest and confounding factors.

_@phy\_tree_ - Phylogenetic tree.

_@refseq_ - A biological sequence set object of a class.


### Example
**Parameter settings:**

```
sample_size <- 100
categories <- 4
OTU_select <- 30
betas0 <- c(1, 0.5, 0.5)
beta <- 1
```

**Generation of microbiome data:**

```
require(phyloseq)
require(multiMiAT)
data("phy", package = "MiHC")
phy.otu.tab <- otu_table(phy)
tree <- phy_tree(phy)
data("phy.parameters")
otu.tab <- SimulateOTU(phy.otu.tab, nSam = 100, phy.parameters, mu = 1000, size = 25)
```

**Random selection of OTUs:**

```
otu.tab_selection <- otu.tab[, order(apply(otu.tab, 2, sum), decreasing = TRUE)[1:OTU_select]]
colnames(otu.tab_selection) <- paste("OTU", colnames(otu.tab_selection), sep = "")
otu.tab_selection <- scale(otu.tab_selection)
```

**Generation of covariates and regression parameters:**

```
x1 <- rbinom(sample_size, 1, 0.5)
x2 <- rnorm(sample_size, 0, 1)
x1 <- scale(x1)
x2 <- scale(x2)
xdata <- data.frame(x1, x2, otu.tab_selection)
beta_coefficients <- c()
for (i in 1:(categories-1)) {
  betas1 <- round(runif(ncol(otu.tab_selection), 0.1, beta), digits = 1)
  beta_coefficients <- c(beta_coefficients, betas0, betas1)
}
beta_coefficients <- c(beta_coefficients, rep(0, (ncol(otu.tab_selection) + 3)))
xnam <- colnames(xdata)
xformula <- as.formula(paste("~ ", paste(xnam, collapse= "+")))
```

**Generation of nomial multicategory outcomes:**

```
y <- nominal.simulated(betas=beta_coefficients, categories=categories, xformula = xformula, xdata = xdata) 
```

**Data synthesis of simulated data with nomial multicategory outcome:**

```
otu.tab_information <- otu_table(otu.tab, taxa_are_rows = FALSE)
sample_information <- sample_data(as.data.frame(cbind(y, x1, x2)))
colnames(sample_information) <- c("y", "x1", "x2")
rownames(sample_information) <- rownames(otu.tab_information)
biom <- phyloseq(otu.tab_information, tree, sample_information)
```


### References
* Liu M, et al. A method for subtype analysis with somatic mutations. _Bioinformatics_ 2021;**37**(1):50-56.

* McMurdie PJ and Holmes S. phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. PLoS ONE. 2013;**8**(4):e61217.

* Touloumis A. Simulating correlated binary and multinomial responses under marginal model specification: The SimCorMultRes package. _The R Journal_ 2016,**8**(2):79.



## Function **sigma_value**

### Description
The median value of distance matrix square (D * D) or distance matrix (D).


### Usage
```
sigma_value(D) 
```


### Arguments
* _D_ - Distance matrix.


### Values
_$sigma.GK_ - Parameter sigma setting of Gaussian kernel function.

_$sigma.EK_ - Parameter sigma setting of exponential kernel function.

_$sigma.LaK_ - Parameter sigma setting of Laplacian kernel function.


### Example
```
require(GUniFrac)
require(multiMiAT)
data(throat.otu.tab)
data(throat.tree)
unifs <- GUniFrac::GUniFrac(throat.otu.tab, throat.tree, alpha = c(0, 0.25, 0.5, 0.75, 1))$unifracs
u.unif <- unifs[, , "d_UW"]
sigma.u.unif.list <- sigma_value(u.unif)
```


### References
* Agarwal D, Zhang NR. Semblance: An empirical similarity kernel on probability spaces. _Science Advances_ 2019;**5**(12):eaau9630.

* Chen J, et al. Associating microbiome composition with environmental covariates using generalized UniFrac distances. _Bioinformatics_ 2012;**28**(16):2106-2113.

* Zhan X, et al. A fast small‐sample kernel independence test for microbiome community‐level association analysis. _Biometrics_ 2017;**73**(4):1453-1463.



## Function **Kernel.Matrix**

### Description
The function, **Kernel.Matrix**, establishes the kernel matrix through single distance and diverse kernel functions (i.e., linear kernel, Gaussian kernel, exponential kernel and Laplacian kernel function).


### Usage
```
Kernel.Matrix(D, kernel, sigma.GK, sigma.EK, sigma.LaK) 
```


### Arguments
* _D_ - Distance matrix. 

* _kernel_ - Kernel function. Default is kernel = c("LK", "GK", "EK", "LaK"). "LK" for linear kernel, "GK" for Gaussian kernel, "EK" for exponential kernel, "LaK" for Laplacian kernel.

* _sigma.GK_ - Parameter setting of Gaussian kernel function.

* _sigma.EK_ - Parameter setting of exponential kernel function.

* _sigma.LaK_ - Parameter setting of Laplacian kernel function.


### Values
_$LK_ - Kernel matrix using linear kernel.

_$GK.sigma1_ - Kernel matrix using Gaussian kernel with sigma = 1.

_$EK.sigma1_ - Kernel matrix using exponential kernel with sigma = 1.

_$LaK.sigma1_ - Kernel matrix using Laplacian kernel with sigma = 1.


### Example
```
require(GUniFrac)
require(multiMiAT)
data(throat.otu.tab)
data(throat.tree)
unifs <- GUniFrac::GUniFrac(throat.otu.tab, throat.tree, alpha = c(0, 0.25, 0.5, 0.75, 1))$unifracs
u.unif <- unifs[, , "d_UW"]
u.unif.k <- Kernel.Matrix(u.unif, kernel = c("LK", "GK", "EK", "LaK"), sigma.GK = 1, sigma.EK = 1, sigma.LaK = 1)
```


### References
* Agarwal D, Zhang NR. Semblance: An empirical similarity kernel on probability spaces. _Science Advances_ 2019;**5**(12):eaau9630.

* Koh H, et al. A distance-based kernel association test based on the generalized linear mixed model for correlated microbiome studies. _Frontiers in Genetics_ 2019;**10**:458.

* Zhan X, et al. A fast small‐sample kernel independence test for microbiome community‐level association analysis. _Biometrics_ 2017;**73**(4):1453-1463.

* Zhao N, et al. Testing in microbiome-profiling studies with MiRKAT, the microbiome regression-based kernel association test. _American Journal of Human Genetics_ 2015;**96**(5):797-807.



## Function **Whole.Kernels**

### Description
The construction of the kernel matrix through diverse distances (i.e., Jaccard dissimilarity, Bray-Curtis distance, Unifrac distance, Generalized UniFrac distance with alpha = 0, Generalized UniFrac distance with alpha = 0.25, Generalized UniFrac distance with alpha = 0.5, Generalized UniFrac distance with alpha = 0.75 and Weighted UniFrac distance) and kernel functions (i.e., linear kernel, Gaussian kernel, exponential kernel and Laplacian kernel function).


### Usage
```
Whole.Kernels(otu.tab, tree, distance = c("Jaccard", "BC", "U.UniFrac", "G.UniFrac.00", "G.UniFrac.25", "G.UniFrac.50", "G.UniFrac.75", "W.UniFrac"), kernel = c("LK", "GK", "EK", "LaK"), sigma.GK = 1, sigma.EK = 1, sigma.LaK = 1)
```


### Arguments
* _otu.tab_ - OTU count table. Rows are samples and columns are OTUs.

* _tree_ - Phylogenetic tree. If NULL, the phylogeny-based distance cannot be calculated.

* _distance_ - Distance. Default is distance = c("Jaccard", "BC", "U.UniFrac", "G.UniFrac.00", "G.UniFrac.25", "G.UniFrac.50", "G.UniFrac.75", "W.UniFrac"). "Jaccard" for Jaccard dissimilarity, "BC" for Bray-Curtis distance, "U.UniFrac" for Unifrac distance, "G.UniFrac.00" for Generalized UniFrac distance with alpha = 0, "G.UniFrac.25" for Generalized UniFrac distance with alpha = 0.25, "G.UniFrac.50" for Generalized UniFrac distance with alpha = 0.5, "G.UniFrac.75" for Generalized UniFrac distance with alpha = 0.75, "W.UniFrac" for Weighted UniFrac distance. In addition, "Jaccard" is an alternative distance.

* _kernel_ - Kernel function. Default is kernel = c("LK", "GK", "EK", "LaK"). "LK" for linear kernel, "GK" for Gaussian kernel, "EK" for exponential kernel, "LaK" for Laplacian kernel. In addition, "EK" is an alternative kernel function.

* _sigma.GK_ - Parameter setting of Gaussian kernel function. Default is sigma.GK = 1.

* _sigma.EK_ - Parameter setting of exponential kernel function. Default is sigma.EK = 1.

* _sigma.LaK_ - Parameter setting of Laplacian kernel function. Default is sigma.LaK = 1.


### Values
_$Jaccard_ - Kernel matrix based on Jaccard dissimilarity and diverse kernel functions (i.e., linear kernel, Gaussian kernel, exponential kernel and Laplacian kernel function).

_$BC_ - Kernel matrix based on Bray-Curtis distance and diverse kernel functions (i.e., linear kernel, Gaussian kernel, exponential kernel and Laplacian kernel function).

_$U.UniFrac_ - Kernel matrix based on Unifrac distance and diverse kernel functions (i.e., linear kernel, Gaussian kernel, exponential kernel and Laplacian kernel function).

_$G.UniFrac.00_ - Kernel matrix based on Generalized UniFrac distance with alpha = 0 and diverse kernel functions (i.e., linear kernel, Gaussian kernel, exponential kernel and Laplacian kernel function).

_$G.UniFrac.25_ - Kernel matrix based on Generalized UniFrac distance with alpha = 0.25 and diverse kernel functions (i.e., linear kernel, Gaussian kernel, exponential kernel and Laplacian kernel function).

_$G.UniFrac.50_ - Kernel matrix based on Generalized UniFrac distance with alpha = 0.5 and diverse kernel functions (i.e., linear kernel, Gaussian kernel, exponential kernel and Laplacian kernel function).

_$G.UniFrac.75_ - Kernel matrix based on Generalized UniFrac distance with alpha = 0.75 and diverse kernel functions (i.e., linear kernel, Gaussian kernel, exponential kernel and Laplacian kernel function).

_$W.UniFrac_ - Kernel matrix based on Weighted UniFrac distance and diverse kernel functions (i.e., linear kernel, Gaussian kernel, exponential kernel and Laplacian kernel function).


### Example
```
require(GUniFrac)
require(ecodist)
require(multiMiAT)
data(throat.otu.tab)
data(throat.tree)
Ks <- Whole.Kernels(throat.otu.tab, throat.tree) 
```


### References

* Bray JR, Curtis JT. An ordination of the upland forest communities of southern Wisconsin. _Ecological Monographs_ 1957;**27**(4):325-349.

* Chen J, et al. Associating microbiome composition with environmental covariates using generalized UniFrac distances. _Bioinformatics_ 2012;**28**(16):2106-2113.

* Jaccard P. The distribution of the flora in the alpine zone. _New Phytologist_ 1912;**11**(2):37-50.

* Koh H, et al. A distance-based kernel association test based on the generalized linear mixed model for correlated microbiome studies. _Frontiers in Genetics_ 2019;**10**:458.

* Lozupone C, et al. Quantitative and qualitative β diversity measures lead to different insights into factors that structure microbial communities. _Applied and Environmental Microbiology_ 2007;**73**(5):1576-1585.

* Lozupone CA, Knight Rob. UniFrac: a new phylogenetic method for comparing microbial communities. _Applied and Environmental Microbiology_ 2007;**73**(5):1576-1585.

* Zhao N, et al. Testing in microbiome-profiling studies with MiRKAT, the microbiome regression-based kernel association test. _American Journal of Human Genetics_ 2015;**96**(5):797-807.



## Function **multiMiRKAT_bcl**

### Description
The function, **multiMiRKAT_bcl**, can be used to detect the association between the microbiome and multicategory outcomes, which is microbiome regression-based kernel association tests based on baseline category logit model.


### Usage
```
multiMiRKAT_bcl(y, Ks, covs = NULL, n.perm = 999, seed = 1)
```


### Arguments
* _y_ - Multicategory outcomes (i.e., host phenotype of interest). Multicategory outcomes mainly contain ordinal multicategory outcomes (e.g., disease severity) nomial multicategory outcomes (e.g., tumor subtype). The data format can be factor, numeric, integer, et al.

* _Ks_ - Kernel matrix based on diverse distances and kernel functions.

* _covs_ - Confounding factors (e.g., age, gender). Default is covs = NULL. The data format must be data frame.

* _n.perm_ - A number of permutations. Default is n.perm = 999. 

* _seed_ - Random number generator. Default is seed = 1.


### Values
_$T.statistic_ - The test statistic of the individual tests based on baseline category logit model.

_$multiMiRKAT.bcl.pvs_ - The _p_-values for the individual tests based on baseline category logit model.

_$multiMiRKAT.bcl.Omnibus.p_ - The _p_-values for the omnibus tests based on baseline category logit model.


### Example
**Import requisite R packages**

```
require(ecodist)
require(GUniFrac)
require(phyloseq)
require(VGAM)
require(MASS)
require(permute)
require(ape)
require(Matrix)
require(multiMiAT)
```

**Import example microbiome data**

```
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
```

**Fit multiMiRKAT_bcl**

```
Ks <- Whole.Kernels(otu.tab, tree) 
out <- multiMiRKAT_bcl(y = y, Ks = Ks, covs = covs)
out
```


### References
* Agarwal D, Zhang NR. Semblance: An empirical similarity kernel on probability spaces. _Science Advances_ 2019;**5**(12):eaau9630.

* Koh H, et al. A powerful microbiome-based association test and a microbial taxa discovery framework for comprehensive association mapping. _Microbiome_ 2017;**5**(1):45.

* Yee TW, et al. The VGAM package for categorical data analysis. _Journal of Statistical Software_ 2010;**32**(10).

* Zhao N, et al. Testing in microbiome-profiling studies with MiRKAT, the microbiome regression-based kernel association test. _American Journal of Human Genetics_ 2015;**96**(5):797-807.



## Function **multiMiRKAT_clm**

### Description
The function, **multiMiRKAT_clm**, can be used to detect the association between the microbiome and multicategory outcomes, which is microbiome regression-based kernel association tests based on cumulative link model.


### Usage
```
multiMiRKAT_clm(y, Ks, covs = NULL, n.perm = 999, seed = 1)
```


### Arguments
* _y_ - Multicategory outcomes (i.e., host phenotype of interest). Multicategory outcomes mainly contain ordinal multicategory outcomes (e.g., disease severity) nomial multicategory outcomes (e.g., tumor subtype). The data format can be factor, numeric, integer, et al.

* _Ks_ - Kernel matrix based on diverse distances and kernel functions.

* _covs_ - Confounding factors (e.g., age, gender). Default is covs = NULL. The data format must be data frame.

* _n.perm_ - A number of permutations. Default is n.perm = 999. 

* _seed_ - Random number generator. Default is seed = 1.


### Values
_$T.statistic_ - The test statistic of the individual tests based on cumulative link model.

_$multiMiRKAT.clm.pvs_ - The _p_-values for the individual tests based on cumulative link model.

_$multiMiRKAT.clm.Omnibus.p_ - The _p_-values for the omnibus tests based on cumulative link model.


### Example
**Import requisite R packages**

```
require(ecodist)
require(GUniFrac)
require(phyloseq)
require(VGAM)
require(MASS)
require(permute)
require(ape)
require(Matrix)
require(multiMiAT)
```

**Import example microbiome data**

```
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
```

**Fit multiMiRKAT_clm**

```
Ks <- Whole.Kernels(otu.tab, tree) 
out <- multiMiRKAT_clm(y = y, Ks = Ks, covs = covs)
out
```


### References
* Agarwal D, Zhang NR. Semblance: An empirical similarity kernel on probability spaces. _Science Advances_ 2019;**5**(12):eaau9630.

* Koh H, et al. A powerful microbiome-based association test and a microbial taxa discovery framework for comprehensive association mapping. _Microbiome_ 2017;**5**(1):45.

* Liu M, et al. A method for subtype analysis with somatic mutations. _Bioinformatics_ 2021;**37**(1):50-56.

* Yee TW, et al. The VGAM package for categorical data analysis. _Journal of Statistical Software_ 2010;**32**(10).

* Zhao N, et al. Testing in microbiome-profiling studies with MiRKAT, the microbiome regression-based kernel association test. _American Journal of Human Genetics_ 2015;**96**(5):797-807.



## Function **multiMiRKAT_crm**

### Description
The function, **multiMiRKAT_crm**, can be used to detect the association between the microbiome and multicategory outcomes, which is microbiome regression-based kernel association tests based on continuation ratio model.


### Usage
```
multiMiRKAT_crm(y, Ks, covs = NULL, n.perm = 999, seed = 1)
```


### Arguments
* _y_ - Multicategory outcomes (i.e., host phenotype of interest). Multicategory outcomes mainly contain ordinal multicategory outcomes (e.g., disease severity) nomial multicategory outcomes (e.g., tumor subtype). The data format can be factor, numeric, integer, et al.

* _Ks_ - Kernel matrix based on diverse distances and kernel functions.

* _covs_ - Confounding factors (e.g., age, gender). Default is covs = NULL. The data format must be data frame.

* _n.perm_ - A number of permutations. Default is n.perm = 999. 

* _seed_ - Random number generator. Default is seed = 1.


### Values
_$T.statistic_ - The test statistic of the individual tests based on continuation ratio model.

_$multiMiRKAT.crm.pvs_ - The _p_-values for the individual tests based on continuation ratio model.

_$multiMiRKAT.crm.Omnibus.p_ - The _p_-values for the omnibus tests based on continuation ratio model.


### Example
**Import requisite R packages**

```
require(ecodist)
require(GUniFrac)
require(phyloseq)
require(VGAM)
require(MASS)
require(permute)
require(ape)
require(Matrix)
require(multiMiAT)
```

**Import example microbiome data**

```
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
```

**Fit multiMiRKAT_crm**

```
Ks <- Whole.Kernels(otu.tab, tree) 
out <- multiMiRKAT_crm(y = y, Ks = Ks, covs = covs)
out
```


### References
* Agarwal D, Zhang NR. Semblance: An empirical similarity kernel on probability spaces. _Science Advances_ 2019;**5**(12):eaau9630.

* Koh H, et al. A powerful microbiome-based association test and a microbial taxa discovery framework for comprehensive association mapping. _Microbiome_ 2017;**5**(1):45.

* Yee TW, et al. The VGAM package for categorical data analysis. _Journal of Statistical Software_ 2010;**32**(10).

* Zhao N, et al. Testing in microbiome-profiling studies with MiRKAT, the microbiome regression-based kernel association test. _American Journal of Human Genetics_ 2015;**96**(5):797-807.



## Function **score.test_bcl**

### Description
The function, **score.test_bcl**, is score test based on baseline category logit model, which is a special case of an individual test in subtype analysis with somatic mutations (SASOM). It can be also used for the association between the microbiome and multicategory outcomes.


### Usage
```
score.test_bcl(y, otu.tab, covs, W) 
```


### Arguments
* _y_ - Multicategory outcomes (i.e., host phenotype of interest). Multicategory outcomes mainly contain ordinal multicategory outcomes (e.g., disease severity) nomial multicategory outcomes (e.g., tumor subtype). The data format can be factor, numeric, integer, et al.

* _otu.tab_ - OTU count table. Rows are samples and columns are OTUs.

* _covs_ - Confounding factors (e.g., age, gender). Default is covs = NULL. The data format must be data frame.

* _W_ - Weight vector of OTUs.


### Values
The _p_-value for the score test based on baseline category logit model.


### Example
**Import requisite R packages:**

```
require(phyloseq)
require(VGAM)
require(MASS)
require(permute)
require(ape)
require(Matrix)
require(multiMiAT)
```

**Import example microbiome data:**

```
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
```

**Fit bcl-based score test:**

```
out <- score.test_bcl(y, otu.tab, covs, W)
```


### References
* Liu M, et al. A method for subtype analysis with somatic mutations. _Bioinformatics_ 2021;**37**(1):50-56.



## Function **MiRKAT_MCN**

### Description
The function, **MiRKAT\_MCN**, is MiRKAT_MC based on cumulative link model, which mainly refers to the **MiRKATMC** function that sets "data.type" to "nominal" in R package, **MiRKATMC**. In addition, we fix a bug in the **MiRKATMC** function, i.e., the p-value combination method HMP reports an error When zero value appears in these combined p-values.


### Usage
```
out.MiRKAT_MCN <- MiRKAT_MCN(y, Ks, covs)
```


### Arguments
* _y_ - Multicategory outcomes (i.e., host phenotype of interest). Multicategory outcomes mainly contain ordinal multicategory outcomes (e.g., disease severity) nomial multicategory outcomes (e.g., tumor subtype). The data format can be factor, numeric, integer, et al.
 
* _Ks_ - Kernel matrix based on diverse distances and kernel functions.

* _covs_ - Confounding factors (e.g., age, gender). Default is covs = NULL. The data format must be data frame.


### Values
_$MiRKAT\_MCN.pvs_ - The _p_-values for the individual tests of MiRKAT-MCN test.

_$MiRKAT\_MCN_ - The _p_-value for the omnibus test of MiRKAT-MCN test.


### Example
**Import requisite R packages:**

```
require(ecodist)
require(GUniFrac)
require(phyloseq)
require(VGAM)
require(MASS)
require(permute)
require(ape)
require(Matrix)
require(multiMiAT)
```

**Import example microbiome data:**

```
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
```

**Fit MiRKAT-MCN:**

```
Ks <- Whole.Kernels(otu.tab, tree, distance = c("BC", "U.UniFrac", "G.UniFrac.00", "G.UniFrac.25", "G.UniFrac.50", "G.UniFrac.75", "W.UniFrac"), kernel = c("LK", "GK", "LaK"))
out.MiRKAT_MCN <- MiRKAT_MCN(y = y, Ks = Ks, covs = covs)
```


### References
* Jiang Z, et al. MiRKAT-MC: A distance-based microbiome kernel association test with multi-categorical outcomes. _Frontiers in Genetics_ 2022;**13**:841764.



## Function **MiRKAT_MCO**

### Description
The function, **MiRKAT\_MCO**, is MiRKAT_MC based on cumulative link model, which mainly refers to the **MiRKATMC** function that sets "data.type" to "ordinal" in R package, **MiRKATMC**. In addition, we fix a bug in the **MiRKATMC** function, i.e., the p-value combination method HMP reports an error When zero value appears in these combined p-values.


### Usage
```
out.MiRKAT_MCO <- MiRKAT_MCO(y, Ks, covs)
```


### Arguments
* _y_ - Multicategory outcomes (i.e., host phenotype of interest). Multicategory outcomes mainly contain ordinal multicategory outcomes (e.g., disease severity) nomial multicategory outcomes (e.g., tumor subtype). The data format can be factor, numeric, integer, et al.

* _Ks_ - Kernel matrix based on diverse distances and kernel functions.

* _covs_ - Confounding factors (e.g., age, gender). Default is covs = NULL. The data format must be data frame.


### Values
_$MiRKAT\_MCO.pvs_ - The _p_-values for the individual tests of MiRKAT-MCO test.

_$MiRKAT\_MCO_ - The _p_-value for the omnibus test of MiRKAT-MCO test.


### Example
**Import requisite R packages:**

```
require(ecodist)
require(GUniFrac)
require(phyloseq)
require(VGAM)
require(MASS)
require(permute)
require(ape)
require(Matrix)
require(multiMiAT)
```

**Import example microbiome data:**

```
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
```

**Fit MiRKAT-MCO:**

```
Ks <- Whole.Kernels(otu.tab, tree, distance = c("BC", "U.UniFrac", "G.UniFrac.00", "G.UniFrac.25", "G.UniFrac.50", "G.UniFrac.75", "W.UniFrac"), kernel = c("LK", "GK", "LaK"))
out.MiRKAT_MCO <- MiRKAT_MCO(y = y, Ks = Ks, covs = covs)
```


### References
* Jiang Z, et al. MiRKAT-MC: A distance-based microbiome kernel association test with multi-categorical outcomes. _Frontiers in Genetics_ 2022;**13**:841764.



## Dataset
biom1: Real data on colorectal cancer.

biom2: Real data on _Clostridium difficile_ infections.

phy.parameters: Estimated parameters of Dirichlet-multinomial distribution.



## Statement
Our code mainly refers to R packages, _MiRKAT_, _OMiAT_, _GLMMMiRKAT_, _SASOM_, _MiSPU_ , _MiRKATMC_ and _GEEMiHC_.
