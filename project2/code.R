#' JIAJUN LI, s2265910
#' Add your own function definitions on this file.

#' neg_log_lik
#
#' @description Evaluate the negated log-likelihood for model A and B
#' @param beta A vector with the beta parameters
#' @param data A `data.frame` with the same variables as the `filament1` data set.
#' Must have columns `CAD_Weight` and `Actual_Weight`
#' @param model Either "A" for a log-linear variance model, or "B" for a proportional
#' scaling error model

neg_log_lik <- function(beta, data, model){
  
  mu <- beta[1] + beta[2]*data[["CAD_Weight"]]
  
  # distinguish between the two models to find the particular standard deviation for the betas
  if(model == "A") {
    sigma <- sqrt(exp(beta[3] + beta[4]*data[["CAD_Weight"]]))
  }else{
    sigma <- sqrt(exp(beta[3])+exp(beta[4]) * (data[["CAD_Weight"]]^2))
  }
  - sum(dnorm(data[["Actual_Weight"]],
              mean = mu,
              sd=sigma,
              log = TRUE))
  
}

#' filament_estimate
#
#' @description Estimate filament models with different variance structure
#' @param data A `data.frame` with the same variables as the `filament1` data set.
#' Must have columns `CAD_Weight` and `Actual_Weight`
#' @param model Either "A" for a log-linear variance model, or "B" for a proportional
#' scaling error model
#' @return An estimation object suitable for use with [filament1_predict()]

filament1_estimate <- function(data, model) {
  model <- match.arg(model, c("A", "B"))
  if (model == "A") {
    beta_start <- c(-0.1, 1.07, -2, 0.05)
  } else {
    beta_start <- c(-0.15, 1.07, -13.5, -6.5)
  }
  opt <- optim(beta_start,
               neg_log_lik,
               data = data,
               model = model,
               hessian = TRUE,
               method = "Nelder-Mead",
               control = list(maxit = 5000)
  )
  fit <- list(
    model = model,
    par = opt$par,
    hessian = opt$hessian
  )
  class(fit) <- c("filament1_estimate", "list")
  fit
}

#' filament1_aux_EV
#' 
#' @description Evaluate the expectation and variance for model A and B
#' @param beta A vector with the beta parameters
#' @param data A `data.frame` containing the required predictors, including `CAD_Weight`
#' @param model Either "A" for a log-linear variance model, or "B" for a proportional
#' scaling error model
#' @param Sigma_beta : If not NULL, an estimate of the covariance matrix for
#                 the uncertainty of estimated betas
#' @return A list with four elements:
#     E : E(y|beta,x)
#     V : Var(y|beta,x)
#     VE : Var(E(y|beta,x)|x) or NULL
#     EV : E(Var(y|beta,x)|x) or NULL

filament1_aux_EV <- function(beta, data, model = c("A", "B"),
                             Sigma_beta = NULL) {
  
  model <- match.arg(model)
  if (model == "A") {
    
    ZE.0 <- model.matrix( ~ 1 + CAD_Weight, data = data)
    ZV.0 <- model.matrix( ~ 1 + CAD_Weight, data = data)
    ZE = cbind(ZE.0, ZV.0 * 0) 
    ZV = cbind(ZE.0 * 0, ZV.0)
    
    VE <- EV <- NULL
    if (!is.null(Sigma_beta)) {
      # E(Var(y|beta,x)|x)
      EV <- exp(ZV %*% beta + rowSums(ZV * (ZV %*% Sigma_beta)) / 2)
      # Var(E(y|beta,x)|x)
      VE <- rowSums(ZE * (ZE %*% Sigma_beta))
    }
    out <- list(
      E = ZE %*% beta,
      V = exp(ZV %*% beta),
      VE = VE,
      EV = EV
    )
  } else {
    
    ZE.0 <- model.matrix( ~ 1 + CAD_Weight, data = data)
    ZV.0 <- model.matrix( ~ 1 + I(CAD_Weight^2), data = data)
    ZE = cbind(ZE.0, ZV.0 * 0) 
    ZV = cbind(ZE.0 * 0, ZV.0)
    
    VE <- EV <- NULL
    if (!is.null(Sigma_beta)) {
      # E(Var(y|beta,x)|x)
      # (pmin: Ignore large Sigma_beta values)
      EV <- ZV %*% exp(beta + pmin(0.5^2, diag(Sigma_beta)) / 2)
      # Var(E(y|beta,x)|x)
      VE <- rowSums(ZE * (ZE %*% Sigma_beta))
    }
    out <- list(
      E = ZE %*% beta,
      V = ZV %*% exp(beta),
      VE = VE,
      EV = EV
    )
  }
  out
}

#' filament1_predict
#
#' @description The summary of predictive distribution for each row of 
#' `newdata` set based on the given `data` set
#' @param data A `data.frame` with the same variables as the `filament1` data set.
#'  Must have columns `CAD_Weight` and `Actual_Weight`.
#' @param newdata the data set used to predict
#' @param model Either "A" for a log-linear variance model, or "B" for a proportional
#' scaling error model
#' @param level significant level for prediction intervals
#' @return a data frame with variables mean,sd,lwr and upr summarizing 
#' predictive distribution for each row of `newdata`.
filament1_predict <- function(data, newdata, model, level){
  # use build in function to estimate the given data
  estimate <- filament1_estimate(data, model)
  # based on the estimate to get the prediction of newdata
  aux_EV <- filament1_aux_EV(estimate$par, newdata, model, solve(estimate$hessian))
  mean <- aux_EV$E
  # get the standard deviation of predictive distribution
  sd <- sqrt(aux_EV$VE + aux_EV$EV)
  lwr <- mean - qt(1 - level/2, df = nrow(data) -2)*sd
  upr <- mean + qt(1 - level/2, df = nrow(data) -2)*sd
  # return the data frame
  tibble(mean, sd, lwr, upr)
}

#' squared_error
#' 
#' @description The function used to compute squared error scores
#' @param pred_dis predictive distribution of the data set with model
#' @param testdata actual value to compute squared error score
#' @return the squared error score
squared_error <- function(pred_dis, testdata){
  (testdata$Actual_Weight - pred_dis$mean)^2
}

#' ds_error
#' 
#' @description The function used to compute Dawid-Sebastiani error scores
#' @param pred_dis predictive distribution of the data set with model
#' @param testdata actual value to compute squared error score
#' @return the Dawid-Sebastiani error score
ds_error <- function(pred_dis, testdata){
  (testdata$Actual_Weight - pred_dis$mean)^2/(pred_dis$sd^2)+log(pred_dis$sd^2)
}

#' leave1out
#' 
#' @description performs leave-one-out cross-validation for the selected model 
#' for each observation
#' @param data A `data.frame` with the same variables as the `filament1` data set.
#' @param model Either "A" for a log-linear variance model, or "B" for a proportional
#' @param level significant level for prediction intervals
#' @return a dataframe with original data and four additional columns 
#' mean, sd, se and ds of leave-one-out prediction means, standard deviations, 
#' and prediction scores.
leave1out <- function(data ,model, level){
  # initialize 4 columns 
  lt_mean <- numeric(nrow(data))
  lt_sd <- numeric(nrow(data))
  score_lt_se <- numeric(nrow(data))
  score_lt_ds <- numeric(nrow(data))
  # loop all rows of the dataframe
  for (i in seq_len(nrow(data))){
    # use the other rows to predict
    pred <- filament1_predict(data = data%>% filter(Index != i), 
                              newdata = data%>%filter(Index == i),
                              model = model , level = level)
    # use the specific row to treat as test data
    obs <- data %>%  filter(Index == i)
    lt_mean[i] <- pred$mean
    lt_sd[i] <- pred$sd
    # calculate two scores
    score_lt_se[i] <- squared_error(pred, obs)
    score_lt_ds[i] <- ds_error(pred, obs)
  }
  data_frame(data, lt_mean, lt_sd, score_lt_se, score_lt_ds)
}

#' arch_loglike
#' 
#' @description evaluates the combined log-likelihood log[p(y|N, φ)] for a collection 
#' y of y-observations
#' @param y a collection y of y-observations
#' @param dataf a data frame with columns `N` and `phi`
#' @return a data frame with the log-likelihood for each row-pair (N, phi)
arch_loglike <- function(y, dataf){
  # get the value of left femur
  y1 <- y[1]
  # get the value of right femur
  y2 <- y[2]
  # initialize the log likelihood column
  arch_log <- numeric(nrow(dataf))
  for (i in seq_len(nrow(dataf))){
    # the combined log-likelihood function
    arch_log[i] <- -lgamma(y1+1) - lgamma(y2+1) -lgamma(dataf$N[i]-y1+1)-
      lgamma(dataf$N[i] - y2 +1) + 2*lgamma(dataf$N[i] +1) +
      (y1 + y2)*log(dataf$phi[i]) + (2*dataf$N[i] - y1 -y2)*log(1-dataf$phi[i])
  }
  data_frame(dataf, arch_log)
}


#' estimate
#' 
#' @description takes inputs y, xi, a, b, and K and implements the Monte Carlo 
#' integration method in question2 to approximate py (y), E(N|y), and E(φ|y).
#' @param y a collection y of y-observations
#' @param xi the prior with Geometric distribution Geom(xi)
#' @param a the prior distribution with Beta(a,b)
#' @param b the prior distribution with Beta(a,b)
#' @param K the number of iteration of Monte Carlo
#' @return the approximate p (y), E(N|y), and E(φ|y) by the Monte Carlo
#' integration method in question2

estimate <- function(y, xi, a, b, K){
  # random set the value of N by its distribution
  N <- rgeom(n = K, prob = xi)
  # random set the value of phi by its distribution
  phi <- rbeta(K, a, b)
  # get the dataframe of row-pair (N, phi)
  data <- data_frame(N, phi)
  # use build function to get the data frame of log-likelihood for each row-pair (N, phi)
  log_p_data <- arch_loglike(y, data)
  # get the log-likelihood for each row-pair (N, phi)
  log_p <- log_p_data$arch_log
  # get the likelihood for each row-pair (N, phi)
  p <- exp(log_p)
  # Monte Carlo estimate p (y)
  p_hat <- sum(p)/K
  # Monte Carlo estimate E(N|y)
  EN <- 1/(p_hat*K) * sum(log_p_data$N * p)
  # Monte Carlo estimate E(φ|y)
  Ephi <- 1/(p_hat*K) * sum(log_p_data$phi * p)
  data_frame(p_hat, EN, Ephi)
}

