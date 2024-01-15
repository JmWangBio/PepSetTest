#' Empirical Bayes moderated t-test
#'
#' Fit a linear model to feature abundance and compute moderated t-statistics via the empirical Bayes method.
#'
#' @param dat a dataframe or matrix of feature (e.g., peptide, protein) abundance
#' @param contrasts.par group levels to be compared separated by dash (e.g., "B-A"
#' if group B is to be compared against group A)
#' @param group list of group levels corresponding to each sample. The order of group
#' levels needs to match that of samples in the feature abundance table.
#' @param logged Boolean variable indicating whether data have been log-transformed
#' @param NPeptide numeric vector indicating number of peptides aggregated for each protein. logNPeptide will be passed to limma-trend. Constant prior variance if null.
#'
#' @return \code{FitContrasts} returns an object of class \code{MArrayLM}. See ?\code{limma::eBayes} for details.
#' @export
#'
#' @author Junmin Wang
#'
#' @references Ritchie, ME, Phipson, B, Wu, D, Hu, Y, Law, CW, Shi, W, and Smyth, GK (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. \emph{Nucleic Acids Research} 43, e47.
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
#' FitContrasts(dat, contrasts.par, group)
#'
FitContrasts <- function(dat, contrasts.par, group, 
                         logged = FALSE, NPeptide = NULL) {
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
  trend <- FALSE
  if (!is.null(NPeptide)) {
    fit$Amean <- log(NPeptide)
    trend <- TRUE
  }
  eBayes.fit <- limma::eBayes(limma::contrasts.fit(fit, contrs),
                              trend = trend)
  return(eBayes.fit)
}
