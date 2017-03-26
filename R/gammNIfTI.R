#' Wrapper to run a Generalized Additive Mixed Effects model on an Nifti and output a parametric map
#'
#' This function is able to run a Generalized Additive Model (GAMM) using the gamm4() function. 
#' The analysis will run in all voxels within the mask and will return parametric and smooth coefficients. 
#' The function will create parametric maps according to the model selected. 
#' The function will return a p-map, t-map, z-map, p-adjusted-map for parametric terms and p-map, z-map, p-adjusted-map for smooth terms.
#' You can select which type of p-value correction you want done on the map 
#' 
#' @param image Input image of type 'nifti' or vector of path(s) to images. If multiple paths, the script will call mergeNifti() and merge across time.
#' @param mask Input mask of type 'nifti' or path to mask. Must be a binary mask
#' @param fourdOut To be passed to mergeNifti, This is the path and file name without the suffix to save the fourd file. Default (NULL) means script won't write out 4D image.
#' @param formula Must be a formula passed to gamm4()
#' @param randomFormula Random effects formual passed to gamm4()
#' @param subjData Dataframe containing all the covariates used for the analysis
#' @param mc.preschedule Argument to be passed to mclapply, whether or not to preschedule the jobs. More info in parallel::mclapply
#' @param ncores Number of cores to use
#' @param method which method of correction for multiple comparisons (default is none)
#' @param outDir Path to the folder where to output parametric maps (Default is Null, only change if you want to write maps out)
#' @param ... Additional arguments passed to gamm4()
#' 
#' @return Returns Parametric maps of the fitted models over the NIfTI image
#' 
#' @export
#' @examples
#' 
#' image <- oro.nifti::nifti(img = array(rnorm(1600, sd=10), dim =c(4,4,4,25)))
#' mask <- oro.nifti::nifti(img = array(c(rep(0,14), rep(1,2)), dim = c(4,4,4,1)))
#' set.seed(1)
#' covs <- data.frame(x = runif(25), y = runif(25), id = rep(1:5,5))
#' fm1 <- "~ s(x) + s(y)"
#' randomFormula <- "~(1|id)"
#' Maps <- gammNIfTI(image, mask, formula = fm1, 
#'                  randomFormula = randomFormula, subjData = covs, ncores = 1,
#'                  method="fdr", REML=TRUE)
#' 
#' 


gammNIfTI <- function (image, mask, fourdOut = NULL, formula, randomFormula, 
                       subjData, mc.preschedule = TRUE, ncores = 1, method = "none", 
                       outDir = NULL, ...) {
  if (missing(image)) {
    stop("image is missing")
  }
  if (missing(mask)) {
    stop("mask is missing")
  }
  if (missing(formula)) {
    stop("formula is missing")
  }
  if (missing(subjData)) {
    stop("subjData is missing")
  }
  if (missing(randomFormula)) {
    stop("randomFormula is missing")
  }
  if (class(formula) != "character") {
    stop("formula class must be character")
  }
  if (class(randomFormula) != "character") {
    stop("randomFormula class must be character")
  }
  models <- vgamm4Param(image, mask, fourdOut = fourdOut, formula = formula, 
                      randomFormula = randomFormula, subjData = subjData, mc.preschedule = mc.preschedule, 
                      ncores = ncores, ...)
  print("Creating parametric maps")
  return(parMap(parameters = models, mask = mask, method = method, 
                outDir = outDir))
}