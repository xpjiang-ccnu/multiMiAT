#' @export
multiMiRKAT_bcl <- function (y, Ks, covs = NULL, n.perm = 999, seed = 1) {
  
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
  n <- length(y.factor)

  if (is.null(covs)) {
    fit <- vglm(y.factor ~ 1, family = multinomial)
    fitprob <- fitted.values(fit)
    prob.vector <- as.vector(t(fitprob[,1:(ncategories - 1)]))
  } else {
    fit <- vglm(y.factor ~ ., data = covs, family = multinomial)
    fitprob <- fitted.values(fit)
    prob.vector <- as.vector(t(fitprob[,1:(ncategories - 1)]))
  }
  
  r <- y.vector - prob.vector
  r.matrix <- matrix(r, n, (ncategories - 1), byrow = T)
  r.s <- r.s.matrix <- list()

  set.seed(seed)
  for (j in 1:n.perm) {
    r.s.matrix[[j]] <- r.matrix[shuffle(n),]
    r.s[[j]] <- matrix(t(r.s.matrix[[j]]), n * (ncategories - 1), 1)
  }
  
  Sigma = cov(r.matrix)
  iSig = solve(Sigma)
  V0 <- kronecker(diag(1, n), iSig)
  
  multiMiRKAT.bcl.pvs <- list()
  allQ.omni <- allQ0.omni <- NULL
  
  if (is.matrix(Ks)){
    Ks <- list(Ks)
  } else if (!is.list(Ks)){
    stop("Please enter either a single kernel matrix or a list of kernel matrices for Ks.")
  }
  
  for (i in 1:length(Ks[[1]])) {
    list.kernels <- list()
    for (ii in 1:length(Ks)) { 
      list.kernels <- c(list.kernels, list(Ks[[ii]][[i]]))
    }
    Qs <- rep(NA, length(list.kernels))
    Q0s <- list()
    for (j in 1:length(list.kernels)) {
      Ks.kronecker <- kronecker(list.kernels[[j]], matrix(1, ncategories - 1, ncategories - 1))
      Qs[j] <- t(r) %*% V0 %*% Ks.kronecker %*% V0 %*% r
      Q0s.inv <- rep(NA, n.perm)
      for (k in 1:n.perm) {
        Q0s.inv[k] <- t(r.s[[k]]) %*% V0 %*% Ks.kronecker %*% V0 %*% r.s[[k]]
      }
      Q0s[[j]] <- Q0s.inv
    }
    mirkat.pvs <- rep(NA, length(list.kernels))
    for (j in 1:length(list.kernels)) {
      mirkat.pvs[j] <- (length(which(abs(Q0s[[j]]) > abs(Qs[[j]]))) + 0.01)/(n.perm + 0.01)
    }
    names(mirkat.pvs) <- names(Ks)
    Q.omni <- min(mirkat.pvs)
    Q0.omni <- rep(NA, n.perm)
    for (l in 1:n.perm) {
      Q0s.n <- list()
      for (m in 1:length(list.kernels)) {
        Q0s.n[[m]] <- Q0s[[m]][-l]
      }
      
      a.Qs <- unlist(lapply(list.kernels, function(x) return(t(r.s[[l]]) %*% V0 %*% kronecker(x, matrix(1, ncategories - 1, ncategories - 1)) %*% V0 %*% r.s[[l]])))
      a.pvs <- unlist(mapply(function(x, y) (length(which(abs(x) > abs(y))) + 0.01)/(n.perm + 0.01), Q0s.n, a.Qs))
      Q0.omni[l] <- min(a.pvs)
    }
    
    p.omni.min <- (length(which(Q0.omni < Q.omni)) + 0.01)/(n.perm + 0.01)
    mirkat.pvs <- c(mirkat.pvs, p.omni.min)
    names(mirkat.pvs)[length(mirkat.pvs)] <- "Omni.MiRKAT"
    allQ0.omni <- cbind(allQ0.omni, Q0.omni)
    allQ.omni <- c(allQ.omni, Q.omni)
    multiMiRKAT.bcl.pvs <- c(multiMiRKAT.bcl.pvs, list(mirkat.pvs))
    
  }
  
  T0 <- apply(allQ0.omni, 1, min)
  T <- min(allQ.omni)
  p.multiMiRKAT.bcl <- (length(which(T0 < T)) + 0.01)/(n.perm + 0.01)
  names(p.multiMiRKAT.bcl) <- "multiMiRKAT-N"
  T.statistic <- list(T = T, T0 = T0)
  names(T.statistic) <- c("T.statistic", "T0.statistic")
  
  multiMiRKAT.bcl.pvs <- lapply(multiMiRKAT.bcl.pvs, function(x) round(x, digits = 4))
  p.multiMiRKAT.bcl <- round(p.multiMiRKAT.bcl, digits = 4)
  names(multiMiRKAT.bcl.pvs) <- names(Ks[[1]])

  multiMiRKAT.bcl <- list(T.statistic = T.statistic, multiMiRKAT.bcl.pvs = multiMiRKAT.bcl.pvs, multiMiRKAT.bcl.Omnibus.p = p.multiMiRKAT.bcl)
  
  return(multiMiRKAT.bcl = multiMiRKAT.bcl)
  
}
