#' build_buckets
#' @description this function makes different buckets for the n_fold cross validation function
#' it takes the train data and put it into differenct buckets. the number of buckets depends on the n_fold number
#'
#' @param train.matrix numeric matrix, training data
#' @param n_folds integer, number of buckets to build
#' @return vector
#' @noRd

build_buckets <- function(train_matrix, n_folds = 5) {
  if (!is.matrix(train_matrix)) {
    stop("train_matrix must be a matrix")
  }

  n_cols <- ncol(train_matrix)

  test_number(
    value = n_folds,
    validation.source = c("build_buckets", "n_folds"),
    min = 2, # minimum 2 folds for cross-validation
    max = floor(n_cols / 3), # maximum folds based on ensuring 3 samples per bucket
    integer_only = TRUE # n_folds must be an integer
  )
  # Check if it's possible to have each bucket with at least 3 samples
  if (n_cols < n_folds * 3) {
    stop("Number of columns in train_matrix is too small to ensure a minimum of 3 samples per bucket.
    You can either increase the sample size of the data or use a smaller number of folds")
  }
  # Calculate the number of columns per fold
  columns_per_fold <- n_cols %/% n_folds
  extra_columns <- n_cols %% n_folds

  # Create the fold indices
  fold_indices <- rep(1:n_folds, each = columns_per_fold)

  # Distribute the extra columns
  if (extra_columns > 0) {
    fold_indices <- c(fold_indices, sample(1:n_folds, extra_columns))
  }

  # Shuffle the fold indices
  index_buckets <- sample(fold_indices)

  # Assign the names
  names(index_buckets) <- colnames(train_matrix)

  return(index_buckets)
}
