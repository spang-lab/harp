#' n-fold cross-validation of the estimated_x algorithm
#'
#' @description performs n-fold cross-validation using the estimated_x algorithm to find the best regularization parameter.
#' @param train_cell_counts matrix representing the cellular compositions of the training samples.
#' @param train_bulks matrix of bulk gene expression for training samples, with genes in rows and samples in columns.
#' @param cell_reference_profile matrix containing the cell reference data, with genes in rows and unique cell types in columns (averaged for each cell type).
#' @param n_folds integer specifying the number of buckets/folds in the cross-validation. for example, if you have 15 samples in your training data and want to
#' perform 5-fold cross-validation, each bucket will have 3 samples.
#' @param iteration single integer specifying the number of iterations.
#' @param lambda_seq sequence of floats representing different values for lambda to find the optimal value. a sequence of values for regularization.
#' @param cv_verbose logical indicating whether to display detailed information during cross-validation.
#' @return list containing the results of the cross-validation.
#' @importFrom DTD estimate_c
#' @noRd


n_fold_cross_validation <- function(
    ...,
    train_cell_counts,
    train_bulks,
    cell_reference_profile,
    n_folds = 5,
    lambda_seq = c(seq(0, 1, by = 0.1), 2^seq(1, 5, by = 1)),
    cv_verbose = TRUE) {
  # safety check for cv_verbose
  check_logical(
    value = cv_verbose,
    validation.source = c("n_fold_cross_validation", "cv_verbose")
  )
  if (cv_verbose) {
    print(
      paste0(
        "Start with cross validation for several lambda at ", Sys.time()
      )
    )
  }

  # safety check for n_folds
  n_cols <- ncol(train_bulks)
  test_number(
    value = n_folds,
    validation.source = c("n_fold_cross_validation", "n_folds"),
    min = 2, # minimum 2 folds for cross-validation
    max = floor(n_cols / 3), # maximum folds based on ensuring 3 samples per bucket
    integer_only = TRUE # n_folds must be an integer
  )
  # safety check for lamda_seq
  validate_lambda_seq(
    lambda_seq = lambda_seq,
    caller_function_name = n_fold_cross_validation,
    allow_single = FALSE
  )

  # making the buckes of cross validation
  set.seed(1)
  bucket_indices <- build_buckets(
    train_matrix = train_bulks,
    n_folds = n_folds
  )

  esti_x_result <- list()
  correlation_results <- list()
  correlation_results_non_negative <- list()
  esti_c_result <- list()

  for (lambda in lambda_seq) {
    if (cv_verbose) {
      pos <- which(lambda_seq == lambda) - 1
      cat(
        "\n calculating for lambda: ", lambda, ", completed ", pos,
        " of ", length(lambda_seq), ", ", 100 * pos / length(lambda_seq), "% \n"
      )
    }
    if (cv_verbose) {
      cat("calculating l_fold: ")
    }
    esti_x_result[[as.character(lambda)]] <- list()
    esti_c_result[[paste0("lambda_", lambda)]] <- list()
    df_cor <- list()
    df_cor_non_negative <- list()

    for (fold in 1:n_folds) {
      if (cv_verbose) {
        cat(fold, "\t")
      }
      # Split the complete training data into test and train:
      test.samples <- names(which(bucket_indices == fold))
      train.samples <- names(which(bucket_indices != fold))

      # reduce the train to only include the cv train samples ...
      train_bulks_reduced <- train_bulks[, train.samples]
      train_cell_counts_reduced <- train_cell_counts[, train.samples]

      # test data set
      test_bulks <- train_bulks[, test.samples]
      test_cell_counts <- train_cell_counts[, test.samples]

      # train a model on the reduced training data
      # estimating the reference profile
      esti_x <- estimate_x(
        ...,
        C = train_cell_counts_reduced,
        Y = train_bulks_reduced,
        X_sc = cell_reference_profile,
        lambda = lambda,
        verbose = TRUE
      )
      esti_x_result[[as.character(lambda)]][[as.character(fold)]] <- esti_x$X
      rm.genes <- rownames(esti_x$X)

      esti_c <- estimate_c(
        X.matrix = esti_x$X,
        new.data = as.matrix(test_bulks[rm.genes, ]),
        DTD.model = rep(1, nrow(esti_x$X)),
        estimate.c.type = "direct"
      )

      # compute correlation for the original C
      df_cor[[as.character(fold)]] <- estimated_c_correlation(
        ...,
        estimated_c = esti_c,
        true_c = test_cell_counts
      )
      if (any(is.na(df_cor[[as.character(fold)]]))) {
        cat("when lambda equals to", lambda, "there is Nas in df_cor. they will be set to zero", "\n")
        df_cor[[as.character(fold)]][is.na(df_cor[[as.character(fold)]])] <- 0
      }
      # store the estimated C result
      esti_c_result[[paste0("lambda_", lambda)]][[as.character(fold)]] <- esti_c
      # this is just putting the negative values in to zeros
      estimated_c_non_negative <- esti_c
      estimated_c_non_negative[estimated_c_non_negative < 0] <- 0

      # correlation for non_negative c
      df_cor_non_negative[[as.character(fold)]] <- estimated_c_correlation(
        ...,
        estimated_c = estimated_c_non_negative,
        true_c = test_cell_counts
      )

      if (any(is.na(df_cor_non_negative[[as.character(fold)]]))) {
        cat("when lambda equals to", lambda, "there is Nas in df_cor_non_negative. they will be set to zero", "\n")
        df_cor_non_negative[[as.character(fold)]][is.na(df_cor_non_negative[[as.character(fold)]])] <- 0
      }
    }
    correlation_results[[paste0("lambda_", lambda)]] <- df_cor
    correlation_results_non_negative[[paste0("lambda_", lambda)]] <- df_cor_non_negative
  }
  correlation_list <- list(
    "correlation" = correlation_results,
    "non_negative_correlation" = correlation_results_non_negative
  )

  # this to get an average correlation (first cell type wise over folds and then over celltypes)
  # for each lambda
  ave_correlation <- lapply(correlation_list, calculate_averages_cv_correlations)

  return_list <- list(
    "correlation" = ave_correlation,
    "estimated_c" = esti_c_result
  )
  return(return_list)
}
