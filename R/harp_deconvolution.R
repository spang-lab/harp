#' harp_deconvolution_model
#'
#' @description this function is a wrapper for the harp_pipeline and returns the main outputs.
#' @param train_data a list containing training bulk expression data with "mixtures" as expression values and the corresponding cell counts as "quantities."
#' @param cell_reference_profile a matrix containing the cell profiles, with genes in rows and unique cell types in columns.
#' @param bulk_data a matrix containing the bulk expression data to be deconvolved, with genes in rows and samples in columns. if this NULL only the reference is estimated
#' @param lambda_seq a numeric value (if no cross-validation is performed) or a sequence of values for regularization.
#' @return a list
#' @export


harp_deconvolution_model <- function(
    ...,
    train_data,
    cell_reference_profile,
    bulk_data = NULL,
    verbose = TRUE) {
    harp_model <- harp_pipeline(
        ...,
        train_data = train_data,
        cell_reference_profile = cell_reference_profile,
        bulk_data = bulk_data,
        n_folds = 5,
        lambda_seq = c(seq(0, 1, by = 0.1), 2^seq(1, 5, by = 1)),
        verbose = TRUE
    )

    return.list <- list(
        "harp_reference_profile" = harp_model$reference_profiles$estimated_reference_second,
        "input_data" = list(
            "initial_cell_reference" = cell_reference_profile,
            "deconvolved_bulk_samples" = bulk_data
        ),
        "harp_model" = harp_model
    )
    if (!is.null(bulk_data)) {
        return.list[["estimated_cell_composition"]] <- harp_model$estimated_c$estimated_c_second
    }
    return.list
}
