#'
#' Self-contained Peptide Set Test Workflow
#'
#' Given peptide abundance and assignment of peptide sequences to proteins,
#' execute the self-contained peptide set test workflow to compute the log2 fold change,
#' p-value, and adjusted p-value of all proteins identified.
#'
#' @inheritParams CompPepSetTestWorkflow
#'
#' @return \code{SelfContPepSetTestWorkflow} returns a dataframe containing the following columns
#' \item{protein}{unique protein identifier}
#' \item{NPeps}{number of peptides}
#' \item{Direction}{direction of change}
#' \item{PValue}{raw p-value}
#' \item{adj.P.Val}{adjusted p-value}
#' \item{logFC}{average log2 fold change of peptides}
#' \item{Up}{number of upregulated peptides}
#' \item{Down}{number of downregulated peptides}
#' @export
#'
#' @author Junmin Wang
#'
#' @references Wu, D, Lim, E, Francois Vaillant, F, Asselin-Labat, M-L, Visvader, JE, and Smyth, GK (2010). ROAST: rotation gene set tests for complex microarray experiments. \emph{Bioinformatics} 26, 2176-2182
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
#' SelfContPepSetTestWorkflow(dat, contrasts.par = contrasts.par,
#' group = group,
#' pep_mapping_tbl = pep_mapping_tbl,
#' logged = FALSE)
#'
SelfContPepSetTestWorkflow <- function(dat, contrasts.par, group, pep_mapping_tbl,
                                       logged = FALSE) {
  if (logged) {
    dat.m <- as.matrix(dat)
  } else {
    dat.m <- as.matrix(log2(dat))
  }
  group <- factor(group)
  design.m <- stats::model.matrix(~ 0 + group)
  colnames(design.m) <- levels(group)
  fit <- limma::lmFit(dat.m, design.m)
  contrs <- limma::makeContrasts(contrasts = contrasts.par,
                                 levels = design.m)
  ## run limma analysis
  eBayes.fit <- limma::eBayes(limma::contrasts.fit(fit, contrs))
  ## convert eBayes.fit to dataframe
  contrasts.res <- EnframeContrastsRes(eBayes.fit = eBayes.fit) %>%
    dplyr::rename("peptide" = "feature")
  ## retrieve peptide indices
  result <- data.frame(peptide = rownames(dat)) %>%
    dplyr::inner_join(contrasts.res, by = "peptide") %>%
    dplyr::inner_join(pep_mapping_tbl,
                      by = "peptide")
  ## extract all proteins
  target_prot_ids <- unique(as.character(result$protein))
  ## obtain peptide indices for each protein
  pep_indices_lst <- lapply(target_prot_ids, function(x) {
    which(result$protein == x)
  })
  ## run self-contained peptide set test
  test_output <- limma::fry(dat.m,
                            pep_indices_lst,
                            design = design.m,
                            contrast = contrs) %>%
    tibble::rownames_to_column(var = "prot_index")
  test_output$prot_index <- as.numeric(sub("set", "", test_output$prot_index))
  test_output$protein <- target_prot_ids[test_output$prot_index]
  test_output <- test_output[order(test_output$prot_index, decreasing = FALSE), ]
  rownames(test_output) <- NULL
  ## calculate average peptide log2 fold change
  logFC_lst <- sapply(pep_indices_lst, function(myIndex)
    mean(result$logFC[myIndex]))
  ## get number of upregulated peptides
  numUp_lst <- sapply(pep_indices_lst, function(myIndex)
    sum(result$logFC[myIndex] > 0))
  ## get number of downregulated peptides
  numDown_lst <- sapply(pep_indices_lst, function(myIndex)
    sum(result$logFC[myIndex] < 0)) 
  test_output <- test_output %>%
    dplyr::mutate(logFC = logFC_lst,
                  Up = numUp_lst,
                  Down = numDown_lst) %>%
    dplyr::select("protein", "NGenes", "Direction", "PValue", "FDR",
                  "logFC", "Up", "Down") %>%
    dplyr::rename("NPeps" = "NGenes",
                  "adj.P.Val" = "FDR")
  return(test_output)
}
