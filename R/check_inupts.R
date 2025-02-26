#' check_input_data
#'
#' saftey check for function in harp
#' This function checks the input data for correctness, ensuring that the matrices for the reference profile (`C`),
#' the bulk expression profile (`Y`), and the scaled expression matrix (`X_sc`) fullfill required structural
#' and naming standards. It also verifies the presence of mutual genes between `Y` and `X_sc`, rearranges the data
#' to align them properly, and confirms that sample and gene names are consistent across all matrices.
#' @param C  matrix representing the reference profile, with row names as cell types and column names as samples.
#' @param Y  matrix representing the bulk expression profile, with row names as genes and column names as samples.
#' @param X_sc  matrix representing the reference profiles of cell types, with row names as genes and column names as cell types.
#'
#' @return list containing the validated and reordered matrices
#' @noRd
check_input_data <- function(C, Y, X_sc) {
    message("Starting validation of the input data...")

    # Check if C is a matrix
    if (!is.matrix(C)) {
        stop("In check_input: C (cell counts) is not a matrix")
    }

    # Check if Y is a matrix
    if (!is.matrix(Y)) {
        stop("In check_input: Y (bulk expression) is not a matrix")
    }
    # Check if Y is a matrix
    if (!is.matrix(X_sc)) {
        stop("In check_input: X_sc (cell reference profile) is not a matrix")
    }

    # Check if column names of Y are in column names of C
    if (!all(colnames(Y) %in% colnames(C)) || !all(colnames(C) %in% colnames(Y))) {
        stop("In check_input: sample names in Y (bulk expression) and C (cell counts)  must match exactly")
    }

    # Check if rownames of C match colnames of X_sc
    if (!all(rownames(C) %in% colnames(X_sc)) || !all(colnames(X_sc) %in% rownames(C))) {
        stop("In check_input: cell type names in C (cell counts)  and X_sc (cell reference profile) must match exactly")
    }
    # Check if rownames of Y are in rownames of X_sc
    if (!any(all(rownames(Y) %in% rownames(X_sc)))) {
        message("Genes in  cell reference profile and bulk expression profile do not completely match; mutual genes will be selected.")
    }

    # Select mutual genes
    gene <- intersect(rownames(X_sc), rownames(Y))

    # Stop if no mutual genes are found
    if (length(gene) == 0) {
        stop("In check_input: No mutual genes found between reference profile and bulk expression profile")
    }

    # Produce a warning if mutual genes are fewer than 500
    if (length(gene) < 500) {
        warning("In check_input: Fewer than 500 mutual genes found between reference profile and bulk profile")
    }

    # reorder samples in Y to match the order in c
    Y <- Y[gene, colnames(C)]
    # reorder cell types of X_sc to match the order in C
    X_sc <- X_sc[gene, rownames(C)]
    message("Input data validation complete.")

    list(C = C, Y = Y, X_sc = X_sc)
}

#' test_number
#' saftey check for function in harp
#'
#' @param value value to be tested for number properties
#' @param validation.source vector of length 2: [1] the calling function's name, [2] the parameter name being tested
#' @param min minimum
#' @param max maximum
#' @param integer_only logical, if TRUE checks for integer, if FALSE allows any non-negative number
#'
#' @return TRUE, if no error is detected, stops with error otherwise
#' @noRd
test_number <- function(value,
                        validation.source,
                        min,
                        max,
                        integer_only = TRUE) {
    error.message <- paste0("In ", validation.source[1], ": ", "'", validation.source[2], "'")

    if (!is.numeric(value) || length(value) != 1) {
        error.message <- paste0(error.message, " is not a single numeric value")
        stop(error.message, call. = FALSE)
    }

    # for non-integer case, ensure value is >= 0
    if (!integer_only && value < 0) {
        error.message <- paste0(error.message, " must be non-negative")
        stop(error.message, call. = FALSE)
    }
    # check for integer if required
    if (integer_only && round(value) != value) {
        error.message <- paste0(error.message, " is not an integer")
        stop(error.message, call. = FALSE)
    }

    if (value < min) {
        error.message <- paste0(error.message, " is below minimal value")
        stop(error.message, call. = FALSE)
    }

    if (value > max) {
        error.message <- paste0(error.message, " is above maximal value")
        stop(error.message, call. = FALSE)
    }

    return(TRUE)
}

#' validate lambda sequence
#'
#' saftey check for function in harp
#' @param lambda_seq single non-negative number or sequence of non-negative numbers (at least two values recommended)
#' @param caller_function_name name of the calling function for error messages
#' @param allow_single logical; if TRUE, lambda_seq can be a single non-negative number.
#'
#'
#' @return TRUE if validation passes, stops with error otherwise
#' @noRd
validate_lambda_seq <- function(lambda_seq, caller_function_name, allow_single = TRUE) {
    # check if caller_function_name is provided
    if (missing(caller_function_name)) {
        stop("caller_function_name must be provided")
    }

    # check if input is numeric vector
    if (!is.numeric(lambda_seq)) {
        stop(paste0("In ", caller_function_name, ": lambda_seq must be numeric"))
    }

    # for single value
    if (length(lambda_seq) == 1) {
        if (!allow_single) {
            stop(paste0("In ", caller_function_name, ": lambda_seq cannot be a single value when allow_single is FALSE"))
        }
        test_number(
            value = lambda_seq,
            validation.source = c(caller_function_name, "lambda_seq"),
            min = 0,
            max = Inf,
            integer_only = FALSE # allow any non-negative number
        )
    }
    # for sequence
    else {
        # check if any values are negative
        if (any(lambda_seq < 0)) {
            stop(paste0("In ", caller_function_name, ": all values in lambda_seq must be non-negative"))
        }

        # check if sequence has at least two values
        if (length(lambda_seq) < 2) {
            stop(paste0("In ", caller_function_name, ": lambda_seq must contain at least two values (more values recommended)"))
        }
    }

    return(TRUE)
}

#' check_logical
#'
#' saftey check for functions in harp
#' @param value value to be tested for number properties
#' @param validation.source vector of length 2: [1] the calling function's name, [2] the parameter name being tested
#'
#' @return TRUE, or it throws an error
#' @noRd
check_logical <- function(value,
                          validation.source) {
    error.message <- paste0("In ", validation.source[1], ": ", "'", validation.source[2], "'")

    if (any(!is.logical(value)) || length(value) != 1) {
        error.message <- paste0(error.message, " must be a single value, either 'TRUE' or 'FALSE' ")
        stop(error.message, call. = FALSE)
    }
    return(TRUE)
}
