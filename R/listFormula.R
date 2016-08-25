#' Create list of Formulas for each voxel 
#'
#' This function is internal. 
#' This function creates list of formulas that will be passed for analysis.
#' @param x  Index of voxels to be analyzed
#' @param formula covariates to be included in the analysis
#' @keywords internal
#' @export
#' @examples
#' 
#' 
#' x <- 1
#' fm1 <- "~ x1"
#' formula <- listFormula(x, formula = fm1)



listFormula <- function(x, formula) {
  
  stats::as.formula(paste(x, formula, sep=""), env = parent.frame(n=3))
  
}