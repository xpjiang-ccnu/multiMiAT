#' @export
Kernel.Matrix <- function(D, kernel, sigma.GK, sigma.EK, sigma.LaK) 
{
  n <- nrow(D)
  sigma.GK.matrix <- matrix(sigma.GK)
  sigma.EK.matrix <- matrix(sigma.EK)
  sigma.LaK.matrix <- matrix(sigma.LaK)
  LK <- GK <- EK <- LaK <- list()

  if ("LK" %in% kernel) {
    centerM <- diag(n) - 1/n
    K.LK <- -0.5 * centerM %*% (D * D) %*% centerM
    eK.LK <- eigen(K.LK, symmetric = TRUE)
    K.LK.reconstruct <- eK.LK$vector %*% diag(pmax(0, eK.LK$values)) %*% t(eK.LK$vector)
    LK <- list(K.LK.reconstruct)
    names(LK) <- "LK"
  } 
  
  if ("GK" %in% kernel){
    GK.matrix <- apply(sigma.GK.matrix, 1, function(x) exp(- D^2 / (2 * x^2)))
    GK.names <- NULL
    for (i in 1:nrow(sigma.GK.matrix)) {
      K.GK <- matrix(GK.matrix[, i], n, n)
      eK.GK <- eigen(K.GK, symmetric = TRUE)
      K.GK.reconstruct <- eK.GK$vector %*% diag(abs(eK.GK$values)) %*% t(eK.GK$vector)
      GK[[i]] <- K.GK.reconstruct
      GK.names <- c(GK.names, paste("GK.sigma", sigma.GK.matrix[i,], sep = ""))
    }
    names(GK) <- GK.names
  }
  
  if ("EK" %in% kernel){
    EK.matrix <- apply(sigma.EK.matrix, 1, function(x) exp(- D / (2 * x^2)))
    EK.names <- NULL
    for (i in 1:nrow(sigma.EK.matrix)) {
      K.EK <- matrix(EK.matrix[, i], n, n)
      eK.EK <- eigen(K.EK, symmetric = TRUE)
      K.EK.reconstruct <- eK.EK$vector %*% diag(abs(eK.EK$values)) %*% t(eK.EK$vector)
      EK[[i]] <- K.EK.reconstruct
      EK.names <- c(EK.names, paste("EK.sigma", sigma.EK.matrix[i,], sep = ""))
    }
    names(EK) <- EK.names
  }
  
  if ("LaK" %in% kernel){
    LaK.matrix <- apply(sigma.LaK.matrix, 1, function(x) exp(- D / x ))
    LaK.names <- NULL
    for (i in 1:nrow(sigma.LaK.matrix)) {
      K.LaK <- matrix(LaK.matrix[, i], n, n)
      eK.LaK <- eigen(K.LaK, symmetric = TRUE)
      K.LaK.reconstruct <- eK.LaK$vector %*% diag(abs(eK.LaK$values)) %*% t(eK.LaK$vector)
      LaK[[i]] <- K.LaK.reconstruct
      LaK.names <- c(LaK.names, paste("LaK.sigma", sigma.LaK.matrix[i,], sep = ""))
    }
    names(LaK) <- LaK.names
  }
  K.list <- c(LK, GK, EK, LaK)
  return(K.list)
}