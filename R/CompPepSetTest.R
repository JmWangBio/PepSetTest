#' Competitive peptide set test
#'
#' Given peptide-wise t-statistics and assignment of peptides to proteins,
#' conduct peptide set tests to assess differential protein expression.
#'
#' @inheritParams TTestwCor
#' @inheritParams AggPeps
#' @param result output from \code{EnframeContrastsRes}
#' @param stat statistics to be used in the peptide set test.
#' Options include "t" (t-statistic) and "logFC" (log2 fold change).
#' @param cor_coef inter-peptide correlation coefficient(s)
#'
#' @return \code{CompPepSetTest} returns a dataframe containing the following columns
#' \item{protein}{unique protein identifier}
#' \item{NPeps}{number of peptides}
#' \item{Correlation}{inter-peptide correlation coefficient}
#' \item{Direction}{direction of change}
#' \item{PValue}{raw p-value}
#' \item{adj.P.Val}{p-value adjusted via the Benjamini-Hochberg method}
#' \item{logFC}{average log2 fold change of peptides}
#' \item{Up}{number of upregulated peptides}
#' \item{Down}{number of downregulated peptides}
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
#' # Generate groups and contrasts
#' group <- c(rep("A", 3), rep("B", 3))
#' contrasts.par <- "B-A"
#'
#' fit.cont <- FitContrasts(dat, contrasts.par, group)
#' cont.res <- EnframeContrastsRes(fit.cont)
#' 
#' # Run peptide set test based on t-statistics and standard deviation
#' CompPepSetTest(cont.res, pep_mapping_tbl, stat = "t",
#' cor_coef = 0, pepC.estim = "sd")
#' 
#' # Run peptide set test based on log2 fold change and median absolute deviation
#' CompPepSetTest(cont.res, pep_mapping_tbl, stat = "logFC",
#' cor_coef = 0, pepC.estim = "mad")
#'
CompPepSetTest <- function(result, pep_mapping_tbl,
                           stat = c("t", "logFC"),
                           cor_coef = 0,
                           pepC.estim = c("sd", "mad")) {
  result <- result |>
    ## remove NA values
    dplyr::filter(!is.na(stat)) |>
    ## merge with peptide mapping table
    dplyr::inner_join(pep_mapping_tbl,
                      by = c("feature" = "peptide")) |>
    dplyr::rename("peptide" = "feature")
  ## extract all proteins
  target_prot_ids <- unique(as.character(result$protein))
  ## obtain peptide indices for each protein
  pep_indices_lst <- lapply(target_prot_ids, function(x) {
    which(result$protein == x)
  })
  names(pep_indices_lst) <- target_prot_ids
  ## calculate average peptide log2 fold change
  logFC_lst <- unlist(lapply(pep_indices_lst, function(myIndex)
    mean(result$logFC[myIndex])))
  ## get number of upregulated peptides
  numUp_lst <- unlist(lapply(pep_indices_lst, function(myIndex)
    sum(result$logFC[myIndex] > 0)))
  ## get number of downregulated peptides
  numDown_lst <- unlist(lapply(pep_indices_lst, function(myIndex)
    sum(result$logFC[myIndex] < 0)))
  ## combine protein, average log2fc, and number of up/downregulated into a dataframe
  reg_tbl <- data.frame(protein = target_prot_ids,
                        logFC = logFC_lst,
                        Up = numUp_lst,
                        Down = numDown_lst)
  if (length(cor_coef) == 1) {
    cor_coef <- rep(cor_coef, length(target_prot_ids))
    names(cor_coef) <- target_prot_ids
  }
  cor_coef_reordered <- cor_coef[target_prot_ids]
  names(cor_coef_reordered) <- target_prot_ids
  # assign zero corr. coef. for proteins with NA values (e.g., proteins with only one peptide)
  cor_coef_reordered[is.na(cor_coef_reordered)] <- 0
  ## t-test
  if (stat == "t") {
    test_output <- TTestwCor(result$t,
                             pep_indices_lst,
                             inter.pep.cor = cor_coef_reordered,
                             pepC.estim = pepC.estim)
  } else {
    test_output <- TTestwCor(result$logFC,
                             pep_indices_lst,
                             inter.pep.cor = cor_coef_reordered,
                             pepC.estim = pepC.estim)
  }
  ## merge with reg_tbl
  test_output <- test_output |>
    tibble::rownames_to_column(var = "protein") |>
    dplyr::inner_join(reg_tbl, by = "protein")
  return(test_output)
}
