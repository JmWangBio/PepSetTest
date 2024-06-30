#' Aggregate peptide abundance values
#'
#' Given peptide abundance and assignment of peptide sequences to proteins,
#' aggregate peptide abundance values into protein abundance values.
#'
#' @param dat a dataframe or matrix of peptide abundance, or a SummarizedExperiment object where 
#' grouping and peptide-protein mapping are provided in colData and rowData, respectively.
#' @param pep_mapping_tbl a table mapping peptides to proteins. Alternatively, it can be the
#' column name of the protein in rowData if dat is a SummarizedExperiment object.
#' @param method method of aggregation. Options including "sum" (summed peptide intensity)
#' and "robreg" (robust regression with M-Estimation).
#' @param logged Boolean variable indicating whether abundance data have been
#' log-transformed
#'
#' @return \code{AggPeps} returns a list containing a matrix of protein abundance values and a vector of number of peptides
#' @export
#'
#' @author Junmin Wang
#'
#' @examples
#' # Generate random peptide data
#' dat <- 2^matrix(rnorm(3000), ncol = 6)
#' colnames(dat) <- paste0("Sample", 1:6)
#' rownames(dat) <- paste0("Peptide", 1:500)
#'
#' # Generate peptide mapping table
#' pep_mapping_tbl <- data.frame(peptide = paste0("Peptide", 1:500),
#' protein = paste0("Protein", rep(1:100, each = 5)))
#'
#' AggPeps(dat, pep_mapping_tbl, method = "sum",
#' logged = FALSE)
#' 
#' # Store data as a SummarizedExperiment object
#' library(tibble)
#' library(SummarizedExperiment)
#' rowData <- pep_mapping_tbl |> column_to_rownames(var = "peptide")
#' dat.nn <- dat
#' rownames(dat.nn) <- NULL
#' colnames(dat.nn) <- NULL
#' dat.se <- SummarizedExperiment(assays = list(int = dat.nn), rowData = rowData)
#'
#' AggPeps(dat.se, pep_mapping_tbl = "protein", method = "sum",
#' logged = FALSE)
#'
AggPeps <- function(dat, pep_mapping_tbl,
                    method = c("sum", "robreg"),
                    logged = c(TRUE, FALSE)) {
  ## extract information from dat if a SummarizedExperiment object
  if (methods::is(dat, "SummarizedExperiment")) {
    pep_mapping_tbl <- data.frame(peptide = rownames(SummarizedExperiment::rowData(dat)), 
                                  protein = SummarizedExperiment::rowData(dat)[[pep_mapping_tbl]])
    dat <- SummarizedExperiment::assay(dat)
  }
  ## convert to original scale
  if (logged) {
    dat.m <- as.matrix(2^dat)
  } else {
    dat.m <- as.matrix(dat)
  }
  ## extract all proteins
  target_prot_ids <- unique(as.character(pep_mapping_tbl$protein))
  ## obtain peptide indices for each protein
  pep_lst <- lapply(target_prot_ids, function(x) {
    pep_mapping_tbl[pep_mapping_tbl$protein == x, "peptide"]
  })
  names(pep_lst) <- target_prot_ids
  ## obtain number of peptides for each protein
  NPeptide_lst <- unlist(lapply(pep_lst, function(x) {
    length(x)
  }))
  names(NPeptide_lst) <- target_prot_ids
  ## aggregate for each protein
  if (method == "sum") {
    prot.dat <- do.call('rbind', lapply(target_prot_ids, function(x) {
      colSums(dat.m[pep_lst[[x]], , drop = FALSE], na.rm = TRUE)
    }))
  } else if (method == "robreg") {
    prot.dat <- do.call('rbind', lapply(target_prot_ids, function(x) {
      RobustReg(dat.m[pep_lst[[x]], , drop = FALSE])
    }))
  }
  rownames(prot.dat) <- target_prot_ids
  if (method != "robreg")
    prot.dat <- log2(prot.dat)
  ## store results to list
  prot.dat.lst <- list()
  prot.dat.lst$int <- prot.dat
  prot.dat.lst$NPeptide <- NPeptide_lst
  return(prot.dat.lst)
}
