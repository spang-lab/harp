#' harp_pipeline

#' @description this script is to train harp to infer all parameteres (lambda, estimated refrence, alpha, g in dtd) and then deconvolve the bulk data
# 1st step, a reference profile is estimated by estimate x with or without cross validation
# 2nd step, training DTD on the bulk and the estimated reference in the 1st step and deconvolve the bulk samples in training set
# 3rd step, using the estimated cell type proportions in 2nd step and input training cellular compositions, alpha is eatimated
# 4th step, the input cellular composition matrix is modified by alpha
# 5th step, estimate x with modified c and train bulk samples
# 6th step, train DTD using new reference and modified c, and deconvolve the bulk samples in the test data(bulk_data)
#' @param train_data a list containing training bulk expression data with "mixtures" as expression values and the corresponding cell counts as "quantities."
#' @param cell_reference_profile a matrix containing the cell profiles, with genes in rows and unique cell types in columns.
#' @param bulk_data a matrix containing the bulk expression data to be deconvolved, with genes in rows and samples in columns.
#' @param lambda_seq a numeric value (if no cross-validation is performed) or a sequence of values for regularization.
#' @param n_folds integer specifying the number of buckets/folds in the cross-validation. for example, if you have 15 samples in your training data and want to
#' perform 5-fold cross-validation, each bucket will have 3 samples.
#' @return return.list list
#' @export
#' @import DTD
harp_pipeline <- function(
    ...,
    train_data,
    cell_reference_profile,
    bulk_data = NULL,
    n_folds = 5,
    lambda_seq = c(seq(0, 1, by = 0.1), 2^seq(1, 5, by = 1)),
    verbose = TRUE) {
    if (verbose) {
        print(paste0("the pipleline starts at", Sys.time()))
    }
    start.time <- Sys.time()
    # estimating the reference profile
    reference.profile.result.first <- estimate_reference_profile(
        ...,
        train_cell_counts = train_data$quantities,
        train_bulks = train_data$mixtures,
        cell_reference_profile = cell_reference_profile,
        n_folds = n_folds,
        lambda_seq = lambda_seq,
        verbose = TRUE
    )

    reference.profile.first <- reference.profile.result.first$matrix_x
    # normalization of the reference
    reference.profile.first <- sweep(
        reference.profile.first,
        2,
        mean(colSums(train_data$mixtures)) / colSums(reference.profile.first),
        FUN = "*"
    )

    # optimal lambda of the first round
    optimal.lambda.first <- reference.profile.result.first$lambda
    # 2nd: train DTD and estimate the proportions
    DTD.x.first <- train_dtd_model(
        reference_profile = reference.profile.first,
        train_data = train_data,
    )

    estimated.c.train <- DTD.x.first$estimated.propotions

    # estimating c of test set. this is just to know after the first run of model and without correction with alpha how the results look like
    if (!is.null(bulk_data)) {
        estimated.c.first <- estimate_cell_proportions(
            reference_profile = reference.profile.first,
            bulk_data = bulk_data,
            dtd_model = DTD.x.first$dtd.model
        )
    }
    # 3th step: estimating the alpha
    alpha <- estimate_alpha(C_true = estimated.c.train, C_lab = train_data$quantities)

    # 4th step: modify the cell counts from experiment
    modified.train.pheno <- alpha * train_data$quantities

    # 5th estimating a new reference
    reference.profile.result.second <- estimate_reference_profile(
        ...,
        train_cell_counts = modified.train.pheno,
        train_bulks = train_data$mixtures,
        cell_reference_profile = cell_reference_profile,
        n_folds = n_folds,
        lambda_seq = lambda_seq,
        verbose = TRUE
    )

    reference.profile.second <- reference.profile.result.second$matrix_x
    optimal.lambda.second <- reference.profile.result.second$lambda


    # rm.genes <- rownames(reference.profile.second)

    training.data <- list(
        "mixtures" = train_data$mixtures,
        "quantities" = modified.train.pheno
    )

    # 6th train DTD and quantify the proportions
    DTD.x.second <- train_dtd_model(
        reference_profile = reference.profile.second,
        train_data = training.data,
    )

    estimated.c.train.second <- DTD.x.second$estimated.propotions

    # estimate the cell counts from the test test
    if (!is.null(bulk_data)) {
        estimated.c.second <- estimate_cell_proportions(
            reference_profile = reference.profile.second,
            bulk_data = bulk_data,
            dtd_model = DTD.x.second$dtd.model
        )
    }
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    time.taken

    # making a list for the output of the function
    # TODO: check if we use the first in our results
    return.list <- list(
        "reference_profiles" = list(
            "estimated_reference_first" = reference.profile.first,
            "estimated_reference_second" = reference.profile.second
        ),
        "cell_reference_profile" = cell_reference_profile,
        "cv" = list(
            "mean_cor_first" = reference.profile.result.first$cv_mean_cor,
            "mean_cor_second" = reference.profile.result.second$cv_mean_cor,
            "lambda_seq" = lambda_seq
        ),
        "lambda" = list(
            "optimal_lambda_first" = optimal.lambda.first,
            "optimal_lambda_second" = optimal.lambda.second
        ),
        "alpha" = alpha,
        "time" = time.taken
    )

    if (!is.null(bulk_data)) {
        return.list[["estimated_c"]] <- list(
            "estimated_c_first" = estimated.c.first,
            "estimated_c_second" = estimated.c.second
        )
    }

    return(return.list)
}
