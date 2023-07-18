#' Enframe result of LIMMA analysis
#'
#' Convert result of LIMMA analysis into a dataframe.
#'
#' @param eBayes.fit output from \code{FitContrasts}. See ?\code{limma::eBayes} for details.
#'
#' @return \code{EnframeContrastsRes} returns a dataframe containing the following columns
#' \item{feature}{unique feature identifier}
#' \item{logFC}{log2 fold change}
#' \item{AveExpr}{average log2 feature expression}
#' \item{t}{t-statistic}
#' \item{P.Value}{raw p-value}
#' \item{adj.P.Val}{p-value adjusted via the Benjamini-Hochberg method}
#' \item{B}{B-statistic}
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
#' # Generate groups and contrasts
#' group <- c(rep("A", 3), rep("B", 3))
#' contrasts.par <- "B-A"
#'
#' fit.cont <- FitContrasts(dat, contrasts.par, group)
#' EnframeContrastsRes(fit.cont)
#'
EnframeContrastsRes <- function(eBayes.fit) {
  result <- limma::topTable(eBayes.fit, number = Inf) %>%
    tibble::rownames_to_column("feature") %>%
    ## remove rows containing NAs
    dplyr::filter(!is.na(t))
  return(result)
}
