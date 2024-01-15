#' Aggregate peptide abundance values
#'
#' Given peptide abundance and assignment of peptide sequences to proteins,
#' aggregate peptide abundance values into protein abundance values.
#'
#' @param dat a dataframe or matrix of peptide abundance
#' @param pep_mapping_tbl a table mapping peptides to proteins. "pep_mapping_tbl" and "dat" should contain the same peptides.
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
AggPeps <- function(dat, pep_mapping_tbl,
                    method = c("sum", "robreg"),
                    logged = c(TRUE, FALSE)) {
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
  NPeptide_lst <- sapply(pep_lst, function(x) {
    length(x)
  })
  names(NPeptide_lst) <- target_prot_ids
  ## aggregate for each protein
  aggfun <- NA
  if (method == "sum") {
    aggfun <- "colSums"
    prot.dat <- t(sapply(target_prot_ids, function(x) {
      eval(parse(text = paste0(aggfun, "(dat.m[pep_lst[[x]], , drop = FALSE], na.rm = TRUE)")))
    }))
  } else if (method == "robreg") {
    aggfun <- "RobustReg"
    prot.dat <- t(sapply(target_prot_ids, function(x) {
      eval(parse(text = paste0(aggfun, "(dat.m[pep_lst[[x]], , drop = FALSE])")))
    }))
  }
  if (method != "robreg")
    prot.dat <- log2(prot.dat)
  ## store results to list
  prot.dat.lst <- list()
  prot.dat.lst$int <- prot.dat
  prot.dat.lst$NPeptide <- NPeptide_lst
  return(prot.dat.lst)
}
