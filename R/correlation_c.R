#' calculate Correlations Between Estimated and True cellular compositions
#' 
#' calculates row-wise correlations between two matrices (estimated and true cell proportions)
#' @param estimated_c a matrix estimated cellular composions
#' @param true_c a matrix of true cellulat compositions with same dimensions as estimated_c
#' @param ... Additional arguments passed to stats::cor()
#' 
#' @return numeric vector containing correlation coefficients
#' @export 
#' @examples
#' # create example matrices
#' estimated <- matrix(rnorm(20), nrow = 4)
#' true <- matrix(rnorm(20), nrow = 4)
#' rownames(estimated) <- paste0("cell", 1:4)
#' rownames(true) <- paste0("cell", 1:4)
#' colnames(estimated) <- paste0("sample", 1:5)
#' colnames(true) <- paste0("sample", 1:5)
#' 
#' # calculate correlations
#' results <- estimated_c_correlation (estimated_c = estimated, true_c = true)
estimated_c_correlation <- function(..., estimated_c, true_c) {
  # Input validation
  if (!is.matrix(estimated_c) || !is.matrix(true_c)) {
    stop("Both estimated_c and true_c must be matrices")
  }
  if (!identical(dim(estimated_c), dim(true_c))) {
    stop("Dimensions of estimated_c and true_c must match")
  }
  
  # check row names
  if (is.null(rownames(estimated_c)) || is.null(rownames(true_c))) {
    stop("Both matrices must have row names")
  }
  if (!identical(rownames(estimated_c), rownames(true_c))) {
    stop("Row names of estimated_c and true_c must match exactly")
  }
  
  # check column names
  if (is.null(colnames(estimated_c)) || is.null(colnames(true_c))) {
    stop("Both matrices must have column names")
  }
  if (!identical(colnames(estimated_c), colnames(true_c))) {
    stop("Column names of estimated_c and true_c must match exactly")
  }
  
  # intiate result vector
  n_rows <- nrow(estimated_c)
  correlation <- numeric(n_rows)
  names(correlation) <- rownames(estimated_c)
  

  tryCatch({
    # calculate standard deviations 
    sd_values <- apply(true_c, 1, stats::sd, na.rm = TRUE)
    
    correlation <- sapply(seq_len(n_rows), function(i) {
      if (sd_values[i] == 0) return(0)
      if (all(is.na(estimated_c[i, ])) || all(is.na(true_c[i, ]))) return(NA)
      stats::cor(estimated_c[i, ], true_c[i, ], use = "complete.obs", ...)
    })
    

    names(correlation) <- rownames(estimated_c)
    
  }, error = function(e) {
    warning("Error in correlation calculation: ", e$message)
    return(rep(NA, n_rows))
  })
  
  return(correlation)
}