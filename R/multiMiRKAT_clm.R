#' @export
multiMiRKAT_clm <- function (y, Ks, covs = NULL, n.perm = 999, seed = 1) {
  
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
  
  y.ordered <- as.ordered(y)
  ncategories <- nlevels(y.ordered)
  y.plus <- rep(y.ordered, each = ncategories - 1)
  Intercept <- rep.int(seq(ncategories - 1), length(y.ordered))
  y.vector <- as.numeric(y.plus == Intercept)
  n <- length(y.ordered)

  if (is.null(covs)) {
    fit <- vglm(y.ordered ~ 1, family = VGAM::cumulative(link = "logitlink", parallel = TRUE))
    fitprob <- fitted.values(fit)
    prob.vector <- as.vector(t(fitprob[,1:(ncategories - 1)]))
  } else {
    fit <- vglm(y.ordered ~ ., family = VGAM::cumulative(link = "logitlink", parallel = TRUE), data = covs)
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
  
  multiMiRKAT.clm.pvs <- list()
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
    multiMiRKAT.clm.pvs <- c(multiMiRKAT.clm.pvs, list(mirkat.pvs))
    
  }
  
  T0 <- apply(allQ0.omni, 1, min)
  T <- min(allQ.omni)
  p.multiMiRKAT.clm <- (length(which(T0 < T)) + 0.01)/(n.perm + 0.01)
  names(p.multiMiRKAT.clm) <- "multiMiRKAT-O"
  T.statistic <- list(T = T, T0 = T0)
  names(T.statistic) <- c("T.statistic", "T0.statistic")
  
  multiMiRKAT.clm.pvs <- lapply(multiMiRKAT.clm.pvs, function(x) round(x, digits = 4))
  p.multiMiRKAT.clm <- round(p.multiMiRKAT.clm, digits = 4)
  names(multiMiRKAT.clm.pvs) <- names(Ks[[1]]) 
  
  multiMiRKAT.clm <- list(T.statistic = T.statistic, multiMiRKAT.clm.pvs = multiMiRKAT.clm.pvs, multiMiRKAT.clm.Omnibus.p = p.multiMiRKAT.clm)

  return(multiMiRKAT.clm = multiMiRKAT.clm)
  
}
