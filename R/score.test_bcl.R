#' @export
score.test_bcl <- function (y, otu.tab, covs, W) {
  
  if (!is.matrix(y) & !is.vector(y) & !is.integer(y) & !is.factor(y)) {
    stop("Error: y must be a matrix, vector, integer or factor!")
  } 
  
  if (is.matrix(y)){
    if(any(rowSums(y)!=1)) {
      stop("Error: summation of each row must equal to 1!")
    }  
    y.inter <- which(y == 1, arr.ind = T)[,2]
  }
  
  if (is.vector(y)){
    if(!is.null(covs) && length(y) != nrow(covs)) {
      stop("Error: The dimension of the response variable is not equal to that of the covariable!")
    } else {
      y.inter <- as.integer(y)
    }
  }
  
  if (is.integer(y)) {
    y.inter <- y
  }
  
  y.factor <- as.factor(y)
  ncategories <- nlevels(y.factor)
  y.plus <- rep(y.factor, each = ncategories - 1)
  Intercept <- rep.int(seq(ncategories - 1), length(y.factor))
  y.vector <- as.numeric(y.plus == Intercept)
  y.matrix <- matrix(y.vector, nrow = length(y.factor), ncol = ncategories - 1, byrow = T)
  y.hat.vector <- as.vector(y.matrix)
  n <- length(y.factor)
  
  if (is.null(covs)) {
    fit <- vglm(y.factor ~ 1, family = multinomial)
    fitprob <- fitted.values(fit)[,1:(ncategories - 1)]
    prob.hat.vector <- as.vector(fitprob)
  } else {
    fit <- vglm(y.factor ~ ., data = covs, family = multinomial)
    fitprob <- fitted.values(fit)[,1:(ncategories - 1)]
    prob.hat.vector <- as.vector(fitprob)
  }
  
  otu.tab <- as.matrix(otu.tab)
  w.otu.tab <- otu.tab %*% t(W)
  if (is.null(covs)) {
    covs <- as.matrix(rep(1, nrow(otu.tab)))
  } else {
    covs <- cbind(rep(1, nrow(covs)), covs)
  }
  covs.matrix <- apply(covs, 2, as.numeric)
  covs.s <- kronecker(diag(ncategories - 1), covs.matrix)
  S <- kronecker(diag(ncategories - 1), t(w.otu.tab))
 
  V <- matrix(0, n * (ncategories - 1), n * (ncategories - 1))
  for (i in 1:(ncategories - 2)) 
    for (j in (i + 1):(ncategories - 1)) 
      diag(V[(n * (i - 1) + 1):(n * i), (n * (j - 1) + 1):(n * j)]) <- -fitprob[, i] * fitprob[, j]
  V <- V + t(V)
  diag(V) <- prob.hat.vector * (1 - prob.hat.vector)

  U.theta <- S %*% (y.hat.vector - prob.hat.vector)
  VX <- V %*% covs.s
  V.theta <- S %*% (V - VX %*% solve(t(covs.s) %*% VX) %*% t(VX)) %*% t(S)
  stat.theta <- t(U.theta) %*% solve(V.theta) %*% (U.theta)
  pval.theta <- pchisq(stat.theta, ncategories - 1, lower.tail = FALSE)

  return(pval.theta)
}

