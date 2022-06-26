#' @export
ordinal.simulated <- function(betas, proportion, xformula = formula(xdata), xdata = parent.frame(), link = "logit") {
 
  if (!is.numeric(proportion)) {
    stop("The type of 'proportion' must is number!")
  } 
  
  categories <- length(proportion)
  if (categories < 3 & min(proportion) < 0) {
    stop("The 'proportion' must be positive, the number of 'proportion' must be greater than or equal to 3, and the type of 'proportion' must is number!")
  }
  
  all.proportion <- sum(proportion)
  cumulative.proportion <- proportion[1]/all.proportion
  linear_predictor_formula <- as.formula(xformula)
  
  if (attr(terms(linear_predictor_formula), "intercept") == 0) {
    linear_predictor_formula <- update(linear_predictor_formula, ~. + 1)
  }
  
  xdata <- data.frame(na.omit(xdata))
  xmat <- model.matrix(linear_predictor_formula, data = xdata)
  xmat <- matrix(xmat[, -1], ncol = ncol(xmat) - 1)
  
  if (length(betas) != ncol(xmat)) {
    stop("The length of 'betas' does not match with the number of covariates")
  }
  
  linear_predictor <- matrix(betas, nrow = nrow(xmat), ncol = ncol(xmat), byrow = TRUE) * xmat
  linear_predictor <- matrix(rowSums(linear_predictor), ncol = 1, byrow = TRUE)
  linear_predictor <- as.matrix(linear_predictor)
  sample_size <- nrow(linear_predictor)
  
  if (length(link) != 1) {
    stop("The length of 'link' must be one")
  }
  
  links <- c("probit", "logit", "cloglog", "cauchit")
  
  if (!is.element(link, links)) {
    stop("'link' must be 'probit','logit','cloglog' and/or 'cauchit'")
  }
  
  distr <- switch(link, probit = "qnorm", logit = "qlogis", cloglog = "qgumbel", cauchit = "qcauchy")
  quantile_functions <- as.character(distr)
  ans <- matrix(rnorm(sample_size), sample_size, 1)
  ans <- pnorm(ans)
  quantile_function <- get(quantile_functions, mode = "function")
  simulated_latent_variables <- quantile_function(ans)
  sample_size <- nrow(linear_predictor)
  u_sim <- simulated_latent_variables - linear_predictor
  simulated_responses <- matrix(0, sample_size, 1)
  
  for (i in 2:(categories - 1)) {
    cumulative.proportion <- c(cumulative.proportion, sum(proportion[1:i])/all.proportion)
  }
  
  u_sim_number <- round(nrow(xdata) * cumulative.proportion)
  intercepts <- u_sim[order(u_sim)[u_sim_number], ]
  intercepts <- matrix(intercepts, 1, categories - 1, TRUE)
  intercepts <- cbind(-Inf, intercepts, Inf)
  simulated_responses <- cut(u_sim, intercepts, labels = FALSE)
  return(y=simulated_responses)
}