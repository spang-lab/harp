#' estimate_alpha
#' @description estimating alpha for representation  of dead cells
#' @param C_lab matrix cellular compositions from experiment,eg facs, flow cytometry,or single cell
#' @param C_true matrix cellular compositions estimated from the deconvoltuion tool, DTD
#' @return alpha_estimate vector
#' @noRd

estimate_alpha <- function(
    C_true,
    C_lab) {
  # validation
  if (!is.matrix(C_true) || !is.matrix(C_lab)) {
    stop("Both C_true and C_lab must be matrices")
  }

  if (!identical(dim(C_true), dim(C_lab))) {
    stop("C_true and C_lab must have identical dimensions")
  }

  if (!identical(rownames(C_true), rownames(C_lab))) {
    stop("C_true and C_lab must have identical rownames")
  }

  if (!identical(colnames(C_true), colnames(C_lab))) {
    stop("C_true and C_lab must have identical colnames")
  }

  if (any(is.na(C_true)) || any(is.na(C_lab))) {
    stop("C_true and C_lab cannot contain NA values")
  }
  # initialize vector
  alpha_estimates <- numeric(nrow(C_true))
  names(alpha_estimates) <- rownames(C_true)
  samples <- colnames(C_true)
  C_lab <- C_lab[rownames(C_true), samples]

  for (i in 1:nrow(C_true)) {
    # Create a data frame with C_lab and C_true for input of lm function
    data_df <- data.frame(
      C_lab = C_lab[i, ],
      C_true = C_true[i, ]
    )

    # where all values are zero
    if (all(data_df$C_lab == 0) || all(data_df$C_true == 0)) {
      alpha_estimates[i] <- 1
      next
    }
    model <- stats::lm(C_true ~ C_lab, data = data_df)

    # extract the estimated coefficient for C_lab
    alpha_estimates[i] <- coef(model)["C_lab"]
  }

  alpha_estimates[alpha_estimates <= 0] <- 1 # set negative and zeros to 1

  return(alpha_estimates)
}
