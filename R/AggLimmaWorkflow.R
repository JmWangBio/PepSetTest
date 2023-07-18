#' Aggregation-based LIMMA workflow
#'
#' Given peptide abundance and assignment of peptide sequences to proteins,
#' execute the aggregation-based LIMMA workflow to compute the log2 fold change,
#' p-value, and adjusted p-value of all proteins identified.
#'
#' @inheritParams AggPeps
#' @param contrasts.par group levels to be compared separated by dash (e.g., "B-A"
#' if group B is to be compared against group A)
#' @param group list of group levels corresponding to each sample. The order of group
#' levels needs to match that of samples in the peptide abundance table.
#'
#' @return \code{AggLimmaWorkflow} returns a dataframe containing the following columns
#' \item{feature}{unique protein identifier}
#' \item{logFC}{log2 fold change}
#' \item{AveExpr}{average log2 protein expression}
#' \item{t}{t-statistic}
#' \item{P.Value}{raw p-value}
#' \item{adj.P.Val}{p-value adjusted via the Benjamini-Hochberg method}
#' \item{B}{B-statistic}
#' @export
#'
#' @author Junmin Wang
#' @references Ritchie, ME, Phipson, B, Wu, D, Hu, Y, Law, CW, Shi, W, and Smyth, GK (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. \emph{Nucleic Acids Research} 43, e47.
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
#' AggLimmaWorkflow(dat, contrasts.par = contrasts.par,
#' group = group,
#' pep_mapping_tbl = pep_mapping_tbl,
#' method = "sum",
#' logged = FALSE)
#'
AggLimmaWorkflow <- function(dat, contrasts.par, group,
                             pep_mapping_tbl,
                             method = c("sum", "robreg"),
                             logged = c(TRUE, FALSE)) {
  ## aggregate peptides
  prot.dat <- AggPeps(dat, pep_mapping_tbl, method, logged = logged)
  ## fit statistical model
  eBayes.fit <- FitContrasts(prot.dat, contrasts.par,
                             group, logged = TRUE)
  ## convert eBayes.fit to dataframe
  contrasts.res <- EnframeContrastsRes(eBayes.fit = eBayes.fit) %>%
    dplyr::rename("protein" = "feature")
  return(contrasts.res)
}
