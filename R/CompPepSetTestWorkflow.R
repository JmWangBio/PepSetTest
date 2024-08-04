#' Competitive Peptide Set Test Workflow
#'
#' Given peptide abundance and assignment of peptide sequences to proteins,
#' execute the competitive peptide set test workflow to compute the log2 fold change,
#' p-value, and adjusted p-value of all proteins identified.
#'
#' @inheritParams AggLimmaWorkflow
#' @inheritParams CompPepSetTest
#' @param dat a dataframe or matrix of peptide abundance (row names should be peptide sequences or peptide IDs), 
#' or a SummarizedExperiment object where grouping and peptide-protein mapping are 
#' provided in colData and rowData, respectively.
#' @param group a vector of group levels corresponding to each sample. Alternatively, it can be the 
#' column name of the group in colData if dat is a SummarizedExperiment object.
#' @param covar covariate matrix. Alternatively, it can be the column names of the covariates
#' in colData if dat is a SummarizedExperiment object.
#' @param pep_mapping_tbl a table mapping peptides to proteins (it should include two columns named 
#' "peptide" and "protein"). Alternatively, it can be the column name of the protein in rowData 
#' if dat is a SummarizedExperiment object.
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
#' \item{adj.P.Val}{p-value adjusted via the Benjamini-Hochberg method}
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
#' stat = "t",
#' correlated = TRUE,
#' equal.correlation = TRUE,
#' pepC.estim = "mad",
#' logged = FALSE)
#'
#' # Store data as a SummarizedExperiment object; add covariates
#' library(tibble)
#' library(SummarizedExperiment)
#' colData <- data.frame(sample = LETTERS[seq_along(group)], group = group, 
#' sex = c("M", "F", "M", "F", "F", "M"), age = 1:6) |> 
#' column_to_rownames(var = "sample")
#' rowData <- pep_mapping_tbl |> column_to_rownames(var = "peptide")
#' dat.nn <- dat
#' rownames(dat.nn) <- NULL
#' colnames(dat.nn) <- NULL
#' dat.se <- SummarizedExperiment(assays = list(int = dat.nn), colData = colData, rowData = rowData)
#'
#' CompPepSetTestWorkflow(dat.se, contrasts.par = contrasts.par,
#' group = "group",
#' pep_mapping_tbl = "protein",
#' covar = c("sex", "age"),
#' stat = "t",
#' correlated = TRUE,
#' equal.correlation = TRUE,
#' pepC.estim = "mad",
#' logged = FALSE)
#'
CompPepSetTestWorkflow <- function(dat, contrasts.par, 
                                   group, 
                                   pep_mapping_tbl,
                                   covar = NULL,
                                   stat = c("t", "logFC"),
                                   correlated = FALSE,
                                   equal.correlation = FALSE,
                                   pepC.estim = c("sd", "mad"),
                                   logged = FALSE) {
  ## extract information from dat if a SummarizedExperiment object
  if (methods::is(dat, "SummarizedExperiment")) {
    group <- SummarizedExperiment::colData(dat)[[group]]
    pep_mapping_tbl <- data.frame(peptide = rownames(SummarizedExperiment::rowData(dat)), 
                                  protein = SummarizedExperiment::rowData(dat)[[pep_mapping_tbl]])
    if (!is.null(covar)) {
      covar <- stats::model.matrix(stats::formula(paste("~", paste(covar, collapse = "+"))),
                                   data = SummarizedExperiment::colData(dat))
      check <- apply(covar, 2, function(x) all(x == 1))
      covar <- as.matrix(covar[, !check])     
    }
    dat <- SummarizedExperiment::assay(dat)
  }
  ## fit statistical model
  eBayes.fit <- FitContrasts(dat, contrasts.par, group, covar = covar, logged = logged)
  ## convert eBayes.fit to dataframe
  contrasts.res <- EnframeContrastsRes(eBayes.fit = eBayes.fit)
  ## construct design matrix
  group <- factor(group)
  design <- stats::model.matrix(~ 0 + group)
  design <- cbind(design, covar)
  ## estimate inter-peptide correlation
  if (correlated) {
    inter.pep.cors <- EstimInterPepCor(dat, design,
                                       pep_mapping_tbl = pep_mapping_tbl,
                                       equal.correlation = equal.correlation,
                                       logged = logged)
  } else {
    inter.pep.cors <- 0
  }
  ## run competitive peptide set test
  test_output <- CompPepSetTest(contrasts.res,
                                pep_mapping_tbl = pep_mapping_tbl,
                                stat = stat,
                                cor_coef = inter.pep.cors,
                                pepC.estim = pepC.estim)
  return(test_output)
}
