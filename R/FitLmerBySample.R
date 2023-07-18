#' Fit a linear mixed model
#'
#' Fit a linear mixed model to the abundance of peptides belonging to one protein
#' and compute the correlation coefficient based on variance components. Sample is
#' treated as a random effect in the mixed model.
#'
#' @param y a matrix of log2-transformed peptide abundance for one protein
#' @param design design matrix
#'
#' @return \code{FitLmerBySample} returns the estimated inter-peptide correlation coefficient.
#' @export
#'
#' @author Junmin Wang and Steven Novick
#'
#' @examples
#' y <- matrix(rnorm(1000*6), 1000, 6)
#' design <- cbind(Intercept = 1, Group = c(0, 0, 0, 1, 1, 1))
#'
#' FitLmerBySample(y, design)
#'
FitLmerBySample <- function(y, design) {
  qrdesign <- qr(design)
  y <- qr.resid(qrdesign, t(y))
  y <- t(y)
  y_long <- reshape2::melt(y)
  colnames(y_long) <- c("feature" , "sample", "val")
  fit <- lme4::lmer(val ~ 1|sample, data = y_long)
  res_dat <- as.data.frame(lme4::VarCorr(fit))
  corr <- res_dat[res_dat$grp == "sample", "vcov"] / (
    res_dat[res_dat$grp == "sample", "vcov"] +
      res_dat[res_dat$grp == "Residual", "vcov"]
  )
  return(corr)
}
