################################################################################
## R script for functions in OptHoldoutSize package                           ##
################################################################################
##
## Sami Haidar-Wehbe, Sam Emerson, James Liley
## October 2021
##

## TODO - new file with functions for learning curve estimation
## TODO - add basic error handlers, add random seeds

################################################################################
## Functions for general simulation of random dataset                         ##
################################################################################

##' Coefficients for imperfect risk score
##'
##' @name gen_dr_coefs
##' @description Generate coefficients corresponding to an imperfect risk score, interpretable
##'  as 'baseline' behaviour in the absence of a risk score
##' @keywords simulation
##' @param coefs Original coefficients
##' @param noise Set to TRUE to add Gaussian noise to coefficients
##' @param num_vars Number of variables at hand for baseline calculation
##' @param max_dr_powers ??
##' @return
##' @examples
##'
##' # Lorem Ipsum
gen_dr_coefs <- function(coefs, noise = TRUE, num_vars = 2, max_dr_powers = 1) {
  if (!noise) {
    return(coefs[1:num_vars])
  }

  coefs_dr <- matrix(nrow = max_dr_powers, ncol = length(coefs))

  if (max_dr_powers - 1){
    for (dr_vars in 1:max_dr_powers)
      # If we are estimating the cost for more than one dr power, we use the ad-hoc
      # formula below for the standard deviation
      coefs_dr[dr_vars, ] <- coefs + rnorm(length(coefs), sd = 2 ** (dr_vars - 2))#1))
  } else {
    coefs_dr[1, ] <- coefs + rnorm(length(coefs), sd = 1)
  }

  return(coefs_dr)
}


##' Generate matrix of random observations
##'
##' @name gen_preds
##' @description Generate matrix of random observations. Observations are unit Gaussian-distributed.
##' @keywords simulation
##' @param nobs Number of observations (samples)
##' @param npreds Number of predictors
##' @param ninters Number of interaction terms, default 0
##' @return Data frame of observations
##' @examples
##'
##' # Lorem Ipsum
gen_preds <- function(nobs, npreds, ninters = 0) {
  # Generate gaussian covariates matrix
  X <- as.data.frame(matrix( rnorm(nobs * npreds), nobs, npreds))

  for (i in 1:ninters) {
    X[, npreds + i] <- X[, i] * X[, i + 1]
  }

  return(X)
}


##' Generate response
##'
##' @name gen_resp
##' @description Generate random outcome (response) according to a ground-truth logistic model
##' @keywords simulation
##' @param X Matrix of observations
##' @param coefs Vector of coefficients for logistic model. If NA, random coefficients are generated. Defaults to NA
##' @param coefs_sd If random coefficients are generated, use this SD (mean 0)
##' @param retprobs If TRUE, return class probability; otherwise, return classes. Defaults to FALSE
##' @return Vector of length dim(X)[1] with outcome classes (if retprobs==FALSE) or outcome probabilities (if retprobs==TRUE)
##' @examples
##'
##' # Lorem Ipsum
gen_resp <- function(X, coefs = NA, coefs_sd = 1, retprobs = FALSE) {
  # For now, combines predictors linearly and applies binomial to find class
  nobs <- dim(X)[1]
  npreds <- dim(X)[2]

  # Generate coefficients for each predictor for the linear combination
  if (any(is.na(coefs))) {
    denom <- npreds ** 0.5
    if (!denom) denom <- 1
    coefs <- rnorm(npreds, sd = coefs_sd/denom)
  }

  # First term: linear combination of X's. Second term: Matrix of noise
  #lin_comb <- rowSums(t(t(X) * coefs) + matrix(rnorm(nobs * npreds), nobs, npreds))
  lin_comb <- rowSums(t(t(X) * coefs))
  probs <- logit(lin_comb)

  if (retprobs) return(probs)
  # Round to get binary classes, assuming cuttoff at 0.5
  #classes <- as.data.frame(round(logit(lin_comb)))
  classes <- rbinom(nobs, 1, probs)
  return(list("classes" = classes, "coefs" = coefs))
}



####### Consider revising - this is not specific to baseline predictions; it could be replaced by a function to just calculate logistic probabilities.
##' Generate responses
##'
##' @name oracle_pred
##' @description Probably for deprecation
##' @keywords simulation
##' @param X Matrix of observations
##' @param coefs Vector of coefficients for logistic model.
##' @param num_vars If noise==FALSE, computes using only first num_vars predictors
##' @param noise If TRUE, uses all predictors
##' @return
##' @examples
##'
##' # Lorem Ipsum
oracle_pred <- function(X, coefs, num_vars = 3, noise = TRUE) {
  # We model Dr behaviour as a logistic regression model that has access to the
  # coefficients, with added noise to them, generating imperfect predictions.
  # There is a double layer of indeterminism: firstly, the coefficients are
  # recalculated with new noise every time the function is called, making each
  # prediction different than previous ones; secondly, the returned classes are
  # sampled from a binomial distribution. This is done this way to simulate
  # the fact that human behaviour is not deterministic.

  # Flag to limit power by adding noise to coefficients or limit access to just a number of them

  if ("Y" %in% colnames(X)) {
    X <- X %>% select(-Y)
  }

  if (!noise) {
    X <- X %>% select(all_of(1:num_vars))
  }

  nobs <- dim(X)[1]

  lin_comb <- rowSums(t(t(X) * coefs))
  probs <- logit(lin_comb)

  return(rbinom(nobs, 1, probs))
}


##' Split data
##'
##' @name split_data
##' @description Split data into holdout and intervention sets
##' @keywords simulation
##' @param X Matrix of observations
##' @param frac Fraction of observations to use for the training set
##' @return Vector of TRUE/FALSE values (randomised) with proportion ``frac`` as TRUE
##' @examples
##'
##' # Lorem Ipsum
split_data <- function(X, frac) {
  # Returns a mask for the observations belonging in the training (holdout) set
  nobs <- dim(X)[1]
  mask <- rep(FALSE, nobs) #sample(FALSE, nobs, replace = TRUE)
  train_lim <- floor(nobs * frac)
  mask[1:train_lim] <- TRUE
  mask <- sample(mask)  # Shuffle to avoid correlation between sizes
  return(mask)
}


##' Train model (wrapper)
##'
##' @name model_train
##' @description Train model using either a GLM or a random forest
##' @keywords simulation
##' @param train_data Data to use for training; assumed to have one binary column called `Y`
##' @param model_family Either 'log_reg' for logistic regression or 'rand_forest' for random forest
##' @param ... Passed to function ``glm()`` or ``ranger()``
##' @return Fitted model of type GLM or Ranger
##' @examples
##'
##' # Lorem Ipsum
model_train <- function(train_data, model_family = "log_reg",...) {
  # Takes training data and the model family, returns a model trained on that data
  # Options:
  # log_reg for logistic regression
  # rand_forest for random forest

  if (model_family == "log_reg"){
    model <- glm(Y ~ ., data = train_data, family = binomial(link = "logit"),...)
  } else if (model_family == "rand_forest") {
    model <- ranger(Y ~ ., data = train_data, probability = TRUE,...)
  }

  return(model)
}


## TODO: update to check class of trained_model rather than model_family variable
##' Make predictions
##'
##' @name model_predict
##' @description Make predictions according to a given model
##' @keywords simulation
##' @param data_test Data for which predictions are to be computed
##' @param trained_model Model for which predictions are to be made
##' @param return_type ??
##' @param threshold ??
##' @param model_family ??
##' @param ... Passed to function ``predict.glm()`` or ``predict.ranger()``
##' @return Vector of predictions
##' @examples
##'
##' # Lorem Ipsum
model_predict <- function(data_test, trained_model, return_type, threshold = NULL, model_family = NULL,...) {
  if (model_family == "log_reg") {
    predictions <- predict(trained_model, newdata = data_test, type = "response",...)
  } else if (model_family == "rand_forest") {
    predictions <- predict(trained_model, data = data_test, type = 'response',...)$predictions[ ,2]
  } else if (is.null(model_family)) {
    stop("model_predict: Please provide a correct model family")
  }

  if (return_type == "class") {
    return(ifelse(predictions > threshold, '1', '0'))
  } else if(return_type == "probs"){
    return(predictions)
  } else stop("model_predict: Wrong return type. Specify in return_type")
}

