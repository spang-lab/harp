#' Compute Loss Function
#'
#' Calculates a loss value based on the squared error between predicted bulks and actual
#' bulk expression (sum((X %*% C - Y)^2)). When lambda > 0, adds a regularization term that
#' penalizes differences between X and X_sc using either softplus function
#' (log1p(exp(diff))) for small differences or absolute values when differences
#' become large (>1e8)
#'
#' @param X matrix  estimated reference from harp
#' @param X_sc matrix intial reference of cell types
#' @param C Coefficient matrix for predictions
#' @param Y matrix bulk expression profiles
#' @param lambda single float Regularization parameter
#'
#' @return single float  total loss
#' @noRd
loss_function <- function(X, X_sc, C, Y, lambda) {
  if (lambda > 0) {
    # if these get too large they produce Inf,
    # so do not use approximation but abs in that case
    # otherwise stick with softplus function
    pos <- log1p(exp(X - X_sc))
    neg <- log1p(exp(X_sc - X))
    if (sum(pos) > 1e8 || sum(neg) > 1e8) {
      loss <- sum((X %*% C - Y)^2) +
        lambda * sum(abs(X - X_sc))
    } else {
      loss <- sum((X %*% C - Y)^2) +
        lambda * sum(pos + neg)
    }
  } else {
    loss <- sum((X %*% C - Y)^2)
  }
  return(loss)
}

#' Calculate Gradient for Loss Function
#'
#' Computes the gradient of the loss function with respect to X. The gradient includes
#' two components: the main gradient term 2 * (X %*% C - Y) %*% t(C) which is the
#' gradient of the squared error term, and when lambda > 0, a regularization gradient
#' term lambda * (1/(1 + exp(X_sc - X)) - 1/(1 + exp(X - X_sc))) that uses sigmoid
#' functions to create a smooth gradient for the difference between X and X_sc.
#'
#' @param X matrix estimated reference from harp
#' @param X_sc matrix  intial reference
#' @param C matrix cellular composition
#' @param Y matrix gene expression data
#' @param lambda  single float regularization parameter. If > 0, adds regularization gradient term
#'
#' @return matrix
#' @noRd
gradient_function <- function(X, X_sc, C, Y, lambda) {
  if (lambda > 0) {
    gradient <- 2 * (X %*% C - Y) %*% t(C) +
      lambda * (1 / (1 + exp(X_sc - X)) - 1 / (1 + exp(X - X_sc)))
  } else {
    gradient <- 2 * (X %*% C - Y) %*% t(C)
  }
  return(gradient)
}


#' estimate_x
#' @description this is the function estimate the estimate x, reference profile from the loss
#' function with given lambda
#' @param C matrix, cellular composition, it's a q*n matrix;
#' q = number of the cell types, n = the number of samples
#' @param Y matrix, bulk gene expression, a g * n matrix; g = the number of the genes
#' @param X_sc matrix, intial reference profiles, g*q matrix
#' @param lambda single float
#' @param iteration single integer
#' @param learning_rate single float
#' @param convergence_threshold single float - this parameter is indicates the level of precision for the loss function decided by dimension of Y
#' #the minimum change required between iterations to consider the model parameters as converged
#' @param adapt_rate single integer this indicates the frequency of adapting the learning rate
#' @param gamma single float numeric hyperparameter (0-1) that
#'  controls the gradient norm threshold for algorithm termination
#' @param verbose logical
#' @return ret_list list
#' @noRd
estimate_x <- function(
    C,
    Y,
    X_sc,
    lambda = 1,
    iteration = 100000,
    learning_rate = 0.01,
    decay_rate = 0.5,
    convergence_threshold = 1e-2,
    adapt_step = 100,
    gamma = 1e-4,
    verbose = TRUE) {
  # start time
  start_time <- Sys.time()

  # safety check for verbose
  check_logical(
    value = verbose,
    validation.source = c("estimate_x", "verbose")
  )
  # saftey check of iteration
  test_number(
    value = iteration,
    validation.source = c("estimate_x", "iteration"),
    min = 100000,
    max = Inf,
    integer_only = TRUE # iteration must be an integer
  )
  # safety check for learning_rate
  test_number(
    value = learning_rate,
    validation.source = c("estimate_x", "learning_rate"),
    min = 0.0001,
    max = 0.9,
    integer_only = FALSE
  )
  # safety check for adapt_step
  test_number(
    value = adapt_step,
    validation.source = c("estimate_x", "adapt_step"),
    min = 50,
    max = iteration, # TODO: a better max?
    integer_only = TRUE # adapt_step must be an integer
  )
  # safety check for gamma
  test_number(
    value = gamma,
    validation.source = c("estimate_x", "gamma"),
    min = 1e-4, # TODO: check these value
    max = 1,
    integer_only = FALSE
  )

  # safety check for convergence_threshold

  test_number(
    value = convergence_threshold,
    validation.source = c("estimate_x", "convergence_threshold"),
    min = 1e-6,
    max = 1e-2,
    integer_only = FALSE
  )

  # set parameters
  convergence_threshold_gradient <- gamma * nrow(Y) * ncol(Y)
  convergence_threshold_loss <- convergence_threshold * nrow(Y) * ncol(Y)

  # safety checks for the cell type, sample and gene names and dim and formate of input data
  inputs <- check_input_data(
    C = C,
    Y = Y,
    X_sc = X_sc
  )
  # checked input data
  C <- inputs$C
  Y <- inputs$Y
  X_sc <- inputs$X_sc

  X <- X_sc
  initial_learning_rate <- learning_rate
  n_loss <- matrix(ncol = 2, nrow = iteration)
  # use armijo line search or direct approach
  armijo_backtracking <- TRUE
  learning_rate <- initial_learning_rate

  message("strating to estimate referece (X) ...")
  # gradient decent
  for (i in 1:iteration) {
    # compute loss and gradient
    current_loss <- loss_function(
      X = X,
      X_sc = X_sc,
      C = C,
      Y = Y,
      lambda = lambda
    )
    current_gradient <- gradient_function(
      X = X,
      X_sc = X_sc,
      C = C,
      Y = Y,
      lambda = lambda
    )

    if (verbose && i %% 1000 == 0) {
      cat("Iteration:", i, "\n")
      cat("Loss:", current_loss, "\n")
      cat("Learning Rate:", learning_rate, "\n")
      cat("Gradient", sum(abs(current_gradient)), "\n")
      cat("-------------------------------", "\n")
    }
    # Option 1: Armijo Backtracking for automated learning rate adaption
    if (armijo_backtracking == TRUE) {
      armijo_steps <- 0
      # after one gradient decent step is finished
      # we start again with initial learning rate
      learning_rate <- initial_learning_rate
      # parameter to turn on the while loop
      armijo <- TRUE
      while (armijo) {
        # compute one step
        X_updated <- X - learning_rate * current_gradient
        loss <- loss_function(
          X = X_updated,
          X_sc = X_sc,
          C = C,
          Y = Y,
          lambda = lambda
        )

        n_loss[i, ] <- c(i, loss)
        # checking if the decrease is large enough
        decent <- loss - current_loss
        threshold <- -1 * learning_rate * gamma * sqrt(sum(current_gradient^2))
        if (decent < threshold) {
          # if decrease is large enough we go one step
          armijo <- FALSE
          X <- X_updated
        } else {
          # if decrease is too small we decrease learning rate
          # and start again from top of the while loop
          learning_rate <- learning_rate * decay_rate
          armijo_steps <- armijo_steps + 1
          if (verbose && armijo_steps %% 1 == 0) {
            cat("Iteration:", i, "\n")
            cat("Loss:", loss, "\n")
            cat("Decent:", decent, "\n")
            cat("Thereshold:", threshold, "\n")
            cat("Learning Rate:", learning_rate, "\n")
            cat("Armijo Steps", armijo_steps, "\n")
            cat("-------------------------------", "\n")
          }
          if (learning_rate < 1e-7) {
            break
          }
        }
      }
      if (learning_rate < 1e-7) {
        # If the learning rate gets too small, we can not decrease the loss
        # further so we stop the minimization here
        print("Breaking because loss cannot be decreased any further.")
        print(paste(
          "Stopping at iteration", i, "with loss", round(current_loss, digits = 2),
          "and gradient", round(sum(abs(current_gradient)), digits = 2)
        ))
        break
      }
      # Option 2: straight forward decrease of learning rate
    } else {
      X <- X - learning_rate * current_gradient

      loss <- sum((X %*% C - Y)^2) +
        lambda * sum((X - X_sc)^2)
      n_loss[i, ] <- c(i, current_loss)
      # adapt learning rate every iteration if loss increases
      if (loss > current_loss) {
        learning_rate <- learning_rate * decay_rate
      }
      # adapt learning rate every adapt_every iterations
      if (i %% adapt_step == 0) {
        learning_rate <- learning_rate * decay_rate
      }
    }
    # check if gradient becomes small
    if (i > 2 && (sum(abs(current_gradient)) < convergence_threshold_gradient)) {
      print(paste(
        "Converged at iteration", i,
        "with loss", round(current_loss, digits = 2),
        "and gradient", round(sum(abs(current_gradient)), digits = 2)
      ))
      break
    }
  }
  if (!any(is.na(n_loss)) && (sum(abs(current_gradient)) > convergence_threshold_gradient)) {
    print(paste0(
      "Gradient did not converge after ", iteration, " iterations. you can increase the iteration ",
      "or use smaller lambda"
    ))
  }
  # end time
  end_time <- Sys.time()

  # Print runtime
  cat("Runtime:", round(difftime(end_time, start_time, units = "secs"), digits = 2), "secs\n")


  # set negative value in X to zero
  X[X < 0] <- 0
  cat("The estimation of the reference profile (x) with lambda =", lambda, "has been completed.\n")
  list(
    "X" = X,
    "n_loss" = n_loss,
    "running_time" = round(difftime(end_time, start_time, units = "secs"), digits = 2)
  )
}
