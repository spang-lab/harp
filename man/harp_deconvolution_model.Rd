% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/harp_deconvolution.R
\name{harp_deconvolution_model}
\alias{harp_deconvolution_model}
\title{harp_deconvolution_model}
\usage{
harp_deconvolution_model(
  ...,
  train_data,
  cell_reference_profile,
  bulk_data = NULL,
  verbose = TRUE
)
}
\arguments{
\item{train_data}{a list containing training bulk expression data with "mixtures" as expression values and the corresponding cell counts as "quantities."}

\item{cell_reference_profile}{a matrix containing the cell profiles, with genes in rows and unique cell types in columns.}

\item{bulk_data}{a matrix containing the bulk expression data to be deconvolved, with genes in rows and samples in columns. if this NULL only the reference is estimated}

\item{lambda_seq}{a numeric value (if no cross-validation is performed) or a sequence of values for regularization.}
}
\value{
a list
}
\description{
this function is a wrapper for the harp_pipeline and returns the main outputs.
}
