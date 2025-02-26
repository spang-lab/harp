#' estimate_reference_profile
#'
#' @description this function performs n-fold cross-validation, selects the optimal lambda value from the cross-validation, and estimates the reference profile.
#' @param train_cell_counts matrix representing the cellular compositions of the training samples.
#' @param train_bulks matrix of bulk gene expression for training samples, with genes in rows and samples in columns.
#' @param cell_reference_profile matrix containing the cell reference data, with genes in rows and unique cell types in columns (averaged for each cell type).
#' @param lambda_seq numeric value (if no cross-validation is performed) or a sequence of values for regularization.
#' @param n_folds integer specifying the number of buckets/folds for cross-validation. for example, if there are 15 samples in the training data and 5-fold cross-validation is performed, each bucket will contain 3 samples.
#' @param verbose logical value indicating whether to display detailed information during execution.
#' @return list containing the reference profile.
#' @noRd

estimate_reference_profile <- function(
    ...,
    train_cell_counts,
    train_bulks,
    cell_reference_profile,
    n_folds = 5,
    lambda_seq = c(seq(0, 1, by = 0.1), 2^seq(1, 5, by = 1)),
    verbose = TRUE) {
    # safety check fot verbose
    check_logical(
        value = verbose,
        validation.source = c("estimate_reference_profile ", "verbose")
    )
    if (verbose) {
        cat("validation of Input data...")
    }
    # safety check for lamda_seq
    validate_lambda_seq(
        lambda_seq = lambda_seq,
        caller_function_name = estimate_reference_profile,
        allow_single = TRUE
    )
    # safety checks for the cell type, sample and gene names and dim and formate of input data
    inputs <- check_input_data(
        C = train_cell_counts,
        Y = train_bulks,
        X_sc = cell_reference_profile
    )

    # checked input data
    train_cell_counts <- inputs$C
    train_bulks <- inputs$Y
    cell_reference_profile <- inputs$X_sc

    # transform lambda_seq into a vector (if it is only a number)
    lambda_seq <- c(lambda_seq)

    # in case that only a scalar is provided, directly apply HARP
    # with that lambda value witohut cross validation
    if (length(lambda_seq) == 1) {
        lambda <- lambda_seq[1]
        if (verbose) {
            cat("the reference profile with regularization parameter,lammbda, equal to,", lambda, "is being calculated")
        }
        estimated.ref.x <- estimate_x(
            ...,
            C = train_cell_counts,
            Y = train_bulks,
            X_sc = cell_reference_profile,
            lambda = lambda,
            iteration = 100000,
            learning_rate = 0.01,
            verbose = verbose
        )

        reference_profile <- list(
            "matrix_x" = estimated.ref.x$X,
            "lambda" = lambda
        )
    } else {
        # in case that lambda_seq is actually a sequence
        # run the cross validation along that sequence of lambda values
        cross_validation_result <- n_fold_cross_validation(
            ...,
            train_cell_counts = train_cell_counts,
            train_bulks = train_bulks,
            cell_reference_profile = cell_reference_profile,
            n_folds = n_folds,
            lambda_seq = lambda_seq,
            cv_verbose = verbose
        )

        message(paste("Performing", n_folds, "-fold cross validation"))

        # selecting the lambda with the highests correlation
        message("The cross-validation summary is:", cross_validation_result$correlation$non_negative_correlation$mean)
        message("This corresponds to lambda:", cross_validation_result$correlation$non_negative_correlation$lambda)

        index.max <- which.max(cross_validation_result$correlation$non_negative_correlation$mean)
        optimal.lambda <- cross_validation_result$correlation$non_negative_correlation$lambda[index.max]

        message("Choosing optimal lambda", optimal.lambda)

        # generating the reference
        estimated.ref.x <- estimate_x(
            ...,
            C = train_cell_counts,
            Y = train_bulks,
            X_sc = cell_reference_profile,
            lambda = optimal.lambda,
            verbose = verbose
        )

        reference_profile <- list(
            "matrix_x" = estimated.ref.x$X,
            "lambda" = optimal.lambda,
            "cv_mean_cor" = cross_validation_result$correlation$non_negative_correlation$mean
        )
    }
    reference_profile
}
