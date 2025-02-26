#' this is function to calculate the average of cross validation
#' @param correlation_data the resutl of cross validation
#'
#' @return data frame
#' @noRd

calculate_averages_cv_correlations <- function(correlation_data) {
    correlation_data <- as.data.frame(correlation_data, check.names = FALSE)
    # removing the extra row
    data <- correlation_data[which(rownames(correlation_data) != "extra"), ]
    # Extracting the part of the column names before the period
    # \\. matches the literal dot
    # ([^.]+) matches one or more characters that are not a dot
    # $ ensures we're matching at the end of the string
    column_names <- gsub("\\.([^.]+)$", "", colnames(data))

    # Creating a new data frame to store the averages
    averages <- data.frame(Cell_Type = rownames(data))

    # Calculating the averages for each group
    for (name in unique(column_names)) {
        cols <- grep(paste0("^", name, "\\."), colnames(data))
        averages[[name]] <- rowMeans(data[, cols, drop = FALSE])
    }
    rownames(averages) <- averages$Cell_Type
    averages$Cell_Type <- NULL
    mean <- apply(averages, 2, function(x) {
        mean(x)
    })
    lambda_value <- as.numeric(gsub(".*_", "", colnames(averages)))
    average <- as.data.frame(t(rbind(averages, "mean" = mean, "lambda" = lambda_value)))
    return(average)
}

#' Train DTD Model
#' @param  reference_profile a matrix
#' @param train_data list containing mixtures and quantities of training samples
#' @param verbose logical
#'
#' @return list
#' @import DTD
#' @noRd

train_dtd_model <- function(reference_profile, train_data, verbose = TRUE, ...) {
    # initialization of the g vector for DTD
    start_tweak <- rep(1, nrow(reference_profile))
    names(start_tweak) <- rownames(reference_profile)

    # Here it is crucial for the reference and bulk profiles to be consistent in terms of gene ordering
    gene.ordering <- rownames(reference_profile)
    train_data <- list(
        "mixtures" = train_data$mixtures[gene.ordering, ],
        "quantities" = train_data$quantities
    )

    dtd.model <- train_deconvolution_model(
        tweak = start_tweak,
        X.matrix = reference_profile,
        train.data.list = train_data,
        test.data.list = NULL,
        estimate.c.type = "direct",
        use.implementation = "cxx",
        lambda.seq = "none",
        cv.verbose = verbose,
        verbose = verbose,
        NORM.FUN = "norm1",
        maxit = 10000,
        ...
    )
    proportions <- estimate_c(
        X.matrix = reference_profile,
        new.data = train_data$mixture,
        DTD.model = dtd.model$best.model$Tweak,
        estimate.c.type = "direct"
    )
    proportions[proportions < 0] <- 0

    list(
        "dtd.model" = dtd.model,
        "estimated.propotions" = proportions
    )
}

#' Estimate Cell Proportions
#' @param  reference_profile a matrix
#' @param  bulk_data matrix to be deconvolute
#' @param dtd_model list output of train_deconvolution_model
#' @param estimate.c.type string to indicated if naive least square,nnls should be used to estimate c
#' @return matrix
#' @importFrom DTD estimate_c
#' @noRd

estimate_cell_proportions <- function(
    reference_profile,
    bulk_data,
    dtd_model,
    estimate.c.type = "direct") {
    # Here it is crucial for the reference and bulk profiles to be consistent in terms of gene ordering
    gene.ordering <- rownames(reference_profile)
    proportions <- estimate_c(
        X.matrix = reference_profile,
        new.data = bulk_data[gene.ordering,],
        DTD.model = dtd_model$best.model$Tweak,
        estimate.c.type = "direct"
    )
    proportions[proportions < 0] <- 0
    proportions
}
