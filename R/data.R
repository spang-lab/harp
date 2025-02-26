#' Sample Deconvolution Data
#'
#' This dataset contains a cell type reference constructed from single cell data of 
#' [The landscape of tumor cell states and ecosystems in diffuse large B cell lymphoma](https://www.cell.com/cancer-cell/fulltext/S1535-6108(21)00451-7)
#' by [Steen et al., 2021].
#' Furthermore it contains bulk data (split on the patient level into train and test) artificially constructed by averaging cell profiles of
#' [Dissecting intratumour heterogeneity of nodal B-cell lymphomas at the transcriptional, genetic and drug-response levels](https://www.nature.com/articles/s41556-020-0532-x)
#' by [Roider et al., 2020]
#' Lastly it contains ground truth cell proportions simulating FACS measurments capturing the abundance of cell types
#' @format A list containing three data frames:
#' \describe{
#'   \item{cell_reference_profile}{The anchor reference that is harmonized by harp}
#'   \item{train_data}{The training data for harp, that means bulk counts and cell proportions}
#'   \item{bulk_counts_test}{Independent test counts for the model to be validated on}
#' }
#' @source Deconvolution data of real world sc datasets to be used for a MWE for harp
"harp_data"