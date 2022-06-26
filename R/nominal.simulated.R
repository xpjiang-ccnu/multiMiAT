#' @export
nominal.simulated <- function(betas, categories, xformula = formula(xdata), xdata = parent.frame()) {
  library(evd)
  linear_predictor_formula <- as.formula(xformula)
  
  if (attr(terms(linear_predictor_formula), "intercept") == 0) {
    linear_predictor_formula <- update(linear_predictor_formula, ~. + 1)
  }
  
  xdata <- data.frame(na.omit(xdata))
  xmat <- model.matrix(linear_predictor_formula, data = xdata)
  xmat <- apply(xmat, 2, function(x) rep(x, each = categories))
  
  if (length(betas) != (categories * ncol(xmat))) {
    stop("The length of 'betas' does not match with the provided covariates")
  }
  
  linear_predictor <- matrix(betas, nrow = nrow(xmat), ncol = ncol(xmat), byrow = TRUE) * xmat
  linear_predictor <- matrix(rowSums(linear_predictor), ncol = categories, byrow = TRUE)
  linear_predictor <- as.matrix(linear_predictor)
  sample_size <- nrow(linear_predictor)
  link <- "cloglog"
  distr <- switch(link, probit = "qnorm", logit = "qlogis", cloglog = "qgumbel", cauchit = "qcauchy")
  quantile_functions <- as.character(distr)
  ans <- matrix(rnorm(sample_size*categories), sample_size, categories)
  ans <- pnorm(ans)
  quantile_function <- get(quantile_functions, mode = "function")
  simulated_latent_variables <- quantile_function(ans)
  sample_size <- nrow(linear_predictor)
  u_sim <- simulated_latent_variables + linear_predictor
  u_sim <- matrix(as.vector(t(u_sim)), sample_size, categories, TRUE)
  simulated_responses <- apply(u_sim, 1, which.max)
  simulated_responses <- matrix(simulated_responses, ncol = 1, byrow = TRUE)
  return(y=simulated_responses)
}