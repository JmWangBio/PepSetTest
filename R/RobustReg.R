#' Robust Regression
#'
#' Estimate protein abundance by fitting a linear model to peptide abundance via robust regression.
#'
#' @inheritParams AggPeps
#' @return \code{RobustReg} returns a numeric vector of estimated protein abundance.
#' @export
#'
#' @author Junmin Wang
#'
#' @references Sticker, A, Goeminne, L, Martens, L, and Clement, L (2020). Robust Summarization and Inference in Proteome-wide Label-free Quantification. \emph{Molecular & Cellular Proteomics} 19, 1209-19.
#' @references Gatto, L, Rainer, J, and Gibb, S (2021). MsCoreUtils: Core Utils for Mass Spectrometry Data. R packageversion 1.4.0. https://github.com/RforMassSpectrometry/MsCoreUtils
#'
#' @examples
#' # Generate random peptide data
#' dat <- 2^matrix(rnorm(600), ncol = 6)
#' colnames(dat) <- paste0("Sample", 1:6)
#' rownames(dat) <- paste0("Peptide", 1:100)
#'
#' RobustReg(dat, logged = FALSE)
#'
RobustReg <- function(dat, logged = FALSE) {
  ## convert to log scale
  if (logged) {
    dat.m <- as.matrix(dat)
  } else {
    dat.m <- as.matrix(log2(dat))
  }
  if (is.null(colnames(dat.m)))
    colnames(dat.m) <- paste0('Sample', 1:ncol(dat.m))
  if (is.null(rownames(dat.m)))
    rownames(dat.m) <- paste0('Feature', 1:nrow(dat.m))
  if (nrow(dat.m) == 1)
    return(dat.m)
  if (nrow(dat.m) >= 200) {
    warning("# of peptides exceeds 200; automatically switching to median aggregation.")
    res <- matrixStats::colMedians(dat.m)
    return(res)
  }
  ## convert to long format
  dat.long <- reshape2::melt(dat.m)
  colnames(dat.long) <- c("feature", "sample", "expression")
  ## remove rows containing NA and infinity
  dat.long <- dat.long[!is.na(dat.long$expression) &
                         !is.infinite(dat.long$expression), ]
  ## convert factor to character
  dat.long$sample <- as.character(dat.long$sample)
  dat.long$feature <- as.character(dat.long$feature)
  res <- tryCatch(
    expr = {
      if (length(unique(dat.long$sample)) == 1) {
        dat.long$sample <- rep(1, length(dat.long$sample))
        fit <- MASS::rlm(expression~-1 + feature,
                         data = dat.long)
      } else {
        fit <- MASS::rlm(expression~-1 + sample + feature,
                         data = dat.long)
      }
      sampleid <- seq_along(unique(dat.long$sample))
      coef <- fit$coefficients[sampleid]
      res <- coef[paste0("sample", colnames(dat.m))]
      names(res) <- colnames(dat.m)
      res
    },
    error = function(e) {
      warning('rlm failed; automatically switching to median aggregation.')
      res <- matrixStats::colMedians(dat.m, na.rm = TRUE)
    }
  )
  return(res)
}
