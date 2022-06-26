#' @export
sigma_value <- function (D) 
{
  D2.hat <- D * D
  D.hat <- D
  sigma.GK <- median(D2.hat)
  sigma.EK <- sigma.LaK <- median(D.hat)
  sigma.list <- list(sigma.GK, sigma.EK, sigma.LaK)
  names(sigma.list) <- c("sigma.GK", "sigma.EK", "sigma.LaK")
  return(sigma.list)
}