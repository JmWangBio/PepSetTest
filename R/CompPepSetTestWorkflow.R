#' Competitive Peptide Set Test Workflow
#'
#' Given peptide abundance and assignment of peptide sequences to proteins,
#' execute the competitive peptide set test workflow to compute the log2 fold change,
#' p-value, and adjusted p-value of all proteins identified.
#'
#' @inheritParams EstimInterPepCor
#' @inheritParams CompPepSetTest
#' @param group list of group levels corresponding to each sample. The order of group
#' levels needs to match that of samples in the feature abundance table.
#' @param correlated Boolean variable indicating whether peptides are assumed to be correlated.
#' If correlated, inter-peptide correlation will be estimated.
#' @param equal.correlation Boolean variable indicating whether all pairwise
#' inter-peptide correlation coefficients are assumed to be equal within a protein.
#' If true, the mixed model approach will be applied; otherwise, the approach
#' described in Wu and Smyth (2012), \emph{Nucleic Acids Research} will be applied. Note
#' that this parameter matters only if "correlated" is set to true.
#'
#' @return \code{CompPepSetTestWorkflow} returns a dataframe containing the following columns
#' \item{protein}{unique protein identifier}
#' \item{NPeps}{number of peptides}
#' \item{Correlation}{inter-peptide correlation coefficient}
#' \item{Direction}{direction of change}
#' \item{PValue}{raw p-value}
#' \item{adj.P.Val}{p-value adjusted via the Benjamini-Yekutieli method}
#' \item{logFC}{average log2 fold change of peptides}
#' \item{Up}{number of upregulated peptides}
#' \item{Down}{number of downregulated peptides}
#' @export
#'
#' @author Junmin Wang
#'
#' @references Wu, D, and Smyth, GK (2012). Camera: a competitive gene set test accounting for inter-gene correlation. \emph{Nucleic Acids Research} 40, e133
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
#' CompPepSetTestWorkflow(dat, contrasts.par = contrasts.par,
#' group = group,
#' pep_mapping_tbl = pep_mapping_tbl,
#' stats = "t",
#' correlated = TRUE,
#' equal.correlation = TRUE,
#' pepC.estim = "mad",
#' logged = FALSE)
#'
CompPepSetTestWorkflow <- function(dat, contrasts.par, group, pep_mapping_tbl,
                                   stats = c("t", "logFC"),
                                   correlated = FALSE,
                                   equal.correlation = FALSE,
                                   pepC.estim = c("sd", "mad"),
                                   logged = FALSE) {
  ## fit statistical model
  eBayes.fit <- FitContrasts(dat, contrasts.par, group, logged = logged)
  ## convert eBayes.fit to dataframe
  contrasts.res <- EnframeContrastsRes(eBayes.fit = eBayes.fit)
  ## estimate inter-peptide correlation
  if (correlated) {
    inter.pep.cors <- EstimInterPepCor(dat, contrasts.par, group,
                                       pep_mapping_tbl = pep_mapping_tbl,
                                       equal.correlation = equal.correlation,
                                       logged = logged)
  } else {
    inter.pep.cors <- 0
  }
  ## run competitive peptide set test
  test_output <- CompPepSetTest(contrasts.res,
                                pep_mapping_tbl = pep_mapping_tbl,
                                stats = stats,
                                cor_coef = inter.pep.cors,
                                pepC.estim = pepC.estim)
  return(test_output)
}
