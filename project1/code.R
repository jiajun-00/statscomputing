#' JIAJUN LI, s2265910
#' Add your own function definitions on this file.

#' Log-Exponential density
#'
#' Compute the density or log-density for a Log-Exponential (LogExp)
#' distribution
#'
#' @param x vector of quantiles
#' @param rate vector of rates
#' @param log logical; if TRUE, the log-density is returned

dlogexp <- function(x, rate = 1, log = FALSE) {
  result <- log(rate) + x - rate * exp(x)
  if (!log) {
    exp(result)
  }
  result
}

#' Log-Sum-Exp
#'
#' Convenience function for computing log(sum(exp(x))) in a
#' numerically stable manner
#'
#' @param x numerical vector

log_sum_exp <- function(x) {
  max_x <- max(x, na.rm = TRUE)
  max_x + log(sum(exp(x - max_x)))
}

#' neg_log_like
#' 
#' The function that calculate the negative log likelihood for
#' specified model
#' 
#' @param beta model parameters
#' @param data a data.frame containing the required variables
#' @param model A or B model
#' 
#' @return the negative log likelihood for specified model
neg_log_like <- function(beta, data, model){
  x <- data$CAD_Weight
  # initialize the negative sum
  negative_sum <- 0
  y <- data$Actual_Weight
  if (model == "A") {
    mean_A <- beta[1] + beta[2] * x
    sd_A <- sqrt(exp(beta[3] + beta[4] * x))
    negative_sum <- -sum(dnorm(x = y, mean = mean_A, sd = sd_A, log = TRUE))
  } else if (model == "B") {
    mean_B <- beta[1] + beta[2] * x
    sd_B <- sqrt(exp(beta[3]) + exp(beta[4] )* (x^2))
    negative_sum <- -sum(dnorm(x = y, mean = mean_B, sd = sd_B, log = TRUE))
  }
  return(negative_sum)
}

#' filament1_estimate
#' 
#' The function that uses the R built in function optim() 
#' and neg_log_like() to estimate the two models A and B 
#' using the filament1 data.
#' 
#' @param data a data.frame with the same variables as the filament1 data set
#' @param model model choice (either A or B)
#' 
#' @return par (the best set of parameters found) and 
#' hessian(the estimate of the Hessian at the solution found)
filament1_estimate <- function(data = filament1, model){
  # initialize the initial values for specified model
  initial <- c(0, 0, 0, 0)
  if (model == "A") {
    # change the initial values for model A
    initial <- c(-0.1, 1.07, -2, 0.05)
  } else if (model == "B") {
    # change the initial values for model B
    initial <- c(-0.15, 1.07, -13.5, -6.5)
  }
  opt <- optim(
    initial,
    neg_log_like,
    data = filament1, model = model,
    method = "BFGS",
    hessian = TRUE
  )
  if (opt$convergence != 0){
    print("Convergence is not zero")
  }
  return(list("par" = opt$par, "hessian" = opt$hessian))
}


#' log_prior_density
#' 
#' The function that evaluate the logarithm of the joint prior density 
#' p(θ) for the four θi parameters.
#' 
#' @param theta theta is the θ parameter vector
#' @param params the vector of γ parameters
#' 
#' @return the logarithm of the joint prior density 
#' p(θ) for the four θi parameters.
log_prior_density <- function(theta, params){
  log_theta1 <- dnorm(theta[1], mean = 0, sd = sqrt(params[1]), log = TRUE)
  log_theta2 <- dnorm(theta[2], mean = 1, sd = sqrt(params[2]), log= TRUE)
  log_theta3 <- dlogexp(theta[3], rate = params[3], log = TRUE)
  log_theta4 <- dlogexp(theta[4], rate = params[4], log = TRUE)
  # the logarithm of the joint prior density of four theta
  joint <- log_theta1 + log_theta2 + log_theta3 + log_theta4
  return(joint)
}

#' log_like
#' 
#' The function that evaluates the observation log-likelihood p(y|θ).
#' 
#' @param theta theta is the θ parameter vector
#' @param x CAD weights
#' @param y Actual weights
#' 
#' @return the observation log-likelihood p(y|θ)
log_like <- function(theta, x, y){
  mean_obs <- theta[1] + theta[2] * x
  sd_obs <- sqrt(exp(theta[3]) + exp(theta[4]) * x^2)
  log_like <- sum(dnorm(x = y, mean = mean_obs, sd = sd_obs, log =TRUE))
  return(log_like)
}

#' log_posterior_density
#' 
#' The function that evaluates the logarithm of the posterior density p(θ|y),
#'  apart from some unevaluated normalisation constant.
#' @param theta theta is the θ parameter vector
#' @param x CAD weights
#' @param y Actual weights
#' @param params the vector of γ parameters
#' 
#' @return the logarithm of the posterior density p(θ|y)
log_posterior_density <- function(theta, x, y, params) {
  log_prior_density(theta, params) + log_like(theta, x, y)
}
 
#' posterior mode
#' 
#' The function that find the mode µ of the log-posterior-density and 
#' evaluates the Hessian at the mode as well as 
#' the inverse of the negated Hessian, S.
#' @param theta_start the start vector in optim() function
#' @param x CAD weights
#' @param y Actual weights
#' @param params the vector of γ parameters
#' 
#' @return the mode µ of the log-posterior-density and 
#' evaluates the Hessian at the mode as well as 
#' the inverse of the negated Hessian, S.
posterior_mode <- function(theta_start, x, y, params){
  opt <- optim(
    theta_start,
    log_posterior_density,
    x = x, y = y, params = params,
    hessian = TRUE,
    # to do maximisation instead of minimisation
    control = list(fnscale = -1)
  )
  if (opt$convergence != 0){
    print("Convergence is not zero")
  }
  return(list(mode = opt$par, hessian = opt$hessian, S = solve(-opt$hessian)))
}

#' do_importance
#' 
#' The function that find beta1, beta2, beta3, beta4, log_weights, 
#' containing the βi samples and normalised log-importance-weights.
#' 
#' @param N the number of samples to generate
#' @param mu the mean vector for the importance distribution
#' @param S the covariance matrix
#' @param x CAD weights
#' @param y Actual weights
#' @param params the vector of γ parameters
#' 
#' @return a data.frame with beta1, beta2, beta3, beta4, log_weights, 
#' containing the βi samples and normalised log-importance-weights.
do_imprtance <- function(N, mu, S, x, y, params){
  # get sample by rmvnorm
  sample <- rmvnorm(N, mu, S)
  # initialize the log posterior list
  log_posterior_list = 1:N
  for (i in 1:N) {
    log_posterior_list[i] <- log_posterior_density(sample[i,], x, y, params)
  }
  # get density of the posterior mode
  pm <- dmvnorm(sample, mean = mu, sigma = S, log = TRUE)
  # get logarithm weights
  log_weights <- log_posterior_list - pm
  # normalize the logarithm weights
  normalized_log_weights <- log_weights - log_sum_exp(log_weights)
  beta1 <- sample[,1]
  beta2 <- sample[,2]
  beta3 <- exp(sample[,3])
  beta4 <- exp(sample[,4])
  data.frame(beta1, beta2, beta3, beta4, normalized_log_weights)
}

#' make_CI
#' 
#' The function that find the weighted quantile of given data
#' 
#' @param x numeric vector whose sample quantiles are wanted
#' @param weights numeric vector of non-negative weights
#' @param probs numeric vector of probabilities with values in [0,1]
#' 
#' @return a 1-row, 2-column data.frame with lower and upper quantile
make_CI <- function(x, weights, probs){
  wq <- wquantile(x, probs = probs, weights = weights)
  data.frame(lower = wq[1], upper = wq[2])
}
#' wquantile 
#'
#' Calculates empirical sample quantiles with optional weights, for given probabilities. 
#' Like in quantile(), the smallest observation corresponds to a probability of 0 and the largest to a probability of 1. 
#' Interpolation between discrete values is done when type=7, as in quantile(). 
#' Use type=1 to only generate quantile values from the raw input samples.
#'
#' @param x numeric vector whose sample quantiles are wanted
#' NA and NaN values are not allowed in numeric vectors unless na.rm is TRUE
#' @param probs numeric vector of probabilities with values in [0,1]
#' @param na.rm logical; if true, any NA and NaN's are removed from x before the quantiles are computed
#' @param type numeric, 1 for no interpolation, or 7, for interpolated quantiles. Default is 7
#' @param weights	 numeric vector of non-negative weights, the same length as x, or NULL. The weights are normalised to sum to 1. If NULL, then wquantile(x) behaves the same as quantile(x), with equal weight for each sample value

wquantile <- function (x, probs = seq(0, 1, 0.25), na.rm = FALSE, type = 7, 
                       weights = NULL, ...) 
{
  if (is.null(weights) || (length(weights) == 1)) {
    weights <- rep(1, length(x))
  }
  stopifnot(all(weights >= 0))
  stopifnot(length(weights) == length(x))
  if (length(x) == 1) {
    return(rep(x, length(probs)))
  }
  n <- length(x)
  q <- numeric(length(probs))
  reorder <- order(x)
  weights <- weights[reorder]
  x <- x[reorder]
  wecdf <- pmin(1, cumsum(weights)/sum(weights))
  if (type == 1) {
  }
  else {
    weights2 <- (weights[-n] + weights[-1])/2
    wecdf2 <- pmin(1, cumsum(weights2)/sum(weights2))
  }
  for (pr_idx in seq_along(probs)) {
    pr <- probs[pr_idx]
    if (pr <= 0) {
      q[pr_idx] <- x[1]
    }
    else if (pr >= 1) {
      q[pr_idx] <- x[n]
    }
    else {
      if (type == 1) {
        j <- 1 + pmax(0, pmin(n - 1, sum(wecdf <= pr)))
        q[pr_idx] <- x[j]
      }
      else {
        j <- 1 + pmax(0, pmin(n - 2, sum(wecdf2 <= pr)))
        g <- (pr - c(0, wecdf2)[j])/(wecdf2[j] - c(0, 
                                                   wecdf2)[j])
        q[pr_idx] <- (1 - g) * x[j] + g * x[j + 1]
      }
    }
  }
  q
}

#' Compute empirical weighted cumulative distribution
#'
#' Version of `ggplot2::stat_ecdf` that adds a `weights` property for each
#' observation, to produce an empirical weighted cumulative distribution function.
#' The empirical cumulative distribution function (ECDF) provides an alternative
#' visualisation of distribution. Compared to other visualisations that rely on
#' density (like [geom_histogram()]), the ECDF doesn't require any
#' tuning parameters and handles both continuous and discrete variables.
#' The downside is that it requires more training to accurately interpret,
#' and the underlying visual tasks are somewhat more challenging.
#'
# @inheritParams layer
# @inheritParams geom_point
#' @param na.rm If `FALSE` (the default), removes missing values with
#'    a warning.  If `TRUE` silently removes missing values.
#' @param n if NULL, do not interpolate. If not NULL, this is the number
#'   of points to interpolate with.
#' @param pad If `TRUE`, pad the ecdf with additional points (-Inf, 0)
#'   and (Inf, 1)
#' @section Computed variables:
#' \describe{
#'   \item{x}{x in data}
#'   \item{y}{cumulative density corresponding x}
#' }
#' @seealso wquantile
#' @export
#' @examples
#' library(ggplot2)
#'
#' n <- 100
#' df <- data.frame(
#'   x = c(rnorm(n, 0, 10), rnorm(n, 0, 10)),
#'   g = gl(2, n),
#'   w = c(rep(1/n, n), sort(runif(n))^sqrt(n))
#' )
#' ggplot(df, aes(x, weights = w)) + stat_ewcdf(geom = "step")
#'
#' # Don't go to positive/negative infinity
#' ggplot(df, aes(x, weights = w)) + stat_ewcdf(geom = "step", pad = FALSE)
#'
#' # Multiple ECDFs
#' ggplot(df, aes(x, colour = g, weights = w)) + stat_ewcdf()
#' ggplot(df, aes(x, colour = g, weights = w)) +
#'   stat_ewcdf() +
#'   facet_wrap(vars(g), ncol = 1)

stat_ewcdf <- function(mapping = NULL, data = NULL,
                       geom = "step", position = "identity",
                       ...,
                       n = NULL,
                       pad = TRUE,
                       na.rm = FALSE,
                       show.legend = NA,
                       inherit.aes = TRUE) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = StatEwcdf,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      n = n,
      pad = pad,
      na.rm = na.rm,
      ...
    )
  )
}


#' @title StatEwcdf ggproto object
#' @name StatEwcdf
#' @rdname StatEwcdf
#' @aliases StatEwcdf
#' @format NULL
#' @usage NULL
#' @export
#' @importFrom ggplot2 aes after_stat has_flipped_aes Stat
NULL

StatEwcdf <- ggplot2::ggproto(
  "StatEwcdf", ggplot2::Stat,
  required_aes = c("x|y", "weights"),
  dropped_aes = c("weights"),     
  
  default_aes = ggplot2::aes(y = ggplot2::after_stat(y)),
  
  setup_params = function(data, params) {
    params$flipped_aes <-
      ggplot2::has_flipped_aes(data,
                               params,
                               main_is_orthogonal = FALSE,
                               main_is_continuous = TRUE)
    
    has_x <- !(is.null(data$x) && is.null(params$x))
    has_y <- !(is.null(data$y) && is.null(params$y))
    if (!has_x && !has_y) {
      rlang::abort("stat_ewcdf() requires an x or y aesthetic.")
    }
    has_weights <- !(is.null(data$weights) && is.null(params$weights))
    #    if (!has_weights) {
    #      rlang::abort("stat_ewcdf() requires a weights aesthetic.")
    #    }
    
    params
  },
  
  compute_group = function(data, scales, n = NULL, pad = TRUE, flipped_aes = FALSE) {
    data <- flip_data(data, flipped_aes)
    # If n is NULL, use raw values; otherwise interpolate
    if (is.null(n)) {
      x <- unique(data$x)
    } else {
      x <- seq(min(data$x), max(data$x), length.out = n)
    }
    
    if (pad) {
      x <- c(-Inf, x, Inf)
    }
    if (is.null(data$weights)) {
      data_ecdf <- ecdf(data$x)(x)
    } else {
      data_ecdf <-
        spatstat.geom::ewcdf(
          data$x,
          weights = data$weights / sum(abs(data$weights)) 
        )(x)
    }
    
    df_ecdf <- vctrs::new_data_frame(list(x = x, y = data_ecdf), n = length(x))
    df_ecdf$flipped_aes <- flipped_aes
    ggplot2::flip_data(df_ecdf, flipped_aes)
  }
)





