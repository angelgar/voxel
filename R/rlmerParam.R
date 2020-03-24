#' Run a Linear Mixed Effects Model on all voxels of a NIfTI image and return parametric coefficients and residuals
#'
#'
#' This function is able to run a Linear Mixed Effect Model using the lmer() function.
#' The analysis will run in all voxels in in the mask and will return the model fit for each voxel.
#' The function relies on lmerTest to create p-values using the Satterthwaite Approximation.
#'
#'
#' @param image Input image of type 'nifti' or vector of path(s) to images. If multiple paths, the script will all mergeNifti() and merge across time.
#' @param mask Input mask of type 'nifti' or path to mask. Must be a binary mask
#' @param fourdOut To be passed to mergeNifti, This is the path and file name without the suffix to save the fourd file. Default (NULL) means script won't write out 4D image.
#' @param formula Must be a formula passed to lmer()
#' @param subjData Dataframe containing all the covariates used for the analysis
#' @param mc.preschedule Argument to be passed to mclapply, whether or not to preschedule the jobs. More info in parallel::mclapply
#' @param ncores Number of cores to use
#' @param ... Additional arguments passed to lmer()
#'
#' @keywords internal
#' @return Return list of parametric and spline coefficients (include standard errors and p-values) fitted to each voxel over the masked images passed to function.
#' @export
#'
#'
#'
#' @examples
#'
#'
#' image <- oro.nifti::nifti(img = array(1:1600, dim =c(4,4,4,25)))
#' mask <- oro.nifti::nifti(img = array(c(rep(0,15),1), dim = c(4,4,4,1)))
#' set.seed(1)
#' covs <- data.frame(x = runif(25), id = rep(1:5,5))
#' fm1 <- "~ x + (1|id)"
#' models <- rlmerParam(image, mask, formula = fm1, subjData = covs, ncores = 1)
#'


rlmerParam <- function(image, mask , fourdOut = NULL, formula, subjData, mc.preschedule = TRUE, ncores = 1, ...) {

  if (missing(image)) { stop("image is missing")}
  if (missing(mask)) { stop("mask is missing")}
  if (missing(formula)) { stop("formula is missing")}
  if (missing(subjData)) { stop("subjData is missing")}

  if (class(formula) != "character") { stop("formula class must be character")}

  if (class(image) == "character" & length(image) == 1) {
    image <- oro.nifti::readNIfTI(fname=image)
  } else if (class(image) == "character" & length(image) > 1) {
    image <- mergeNiftis(inputPaths = image, direction = "t", outfile = fourdOut)
  }

  if (class(mask) == "character" & length(mask) == 1) {
    mask <- oro.nifti::readNIfTI(fname=mask)
  }


  imageMat <- ts2matrix(image, mask)


  voxNames <- as.character(names(imageMat))

  rm(image)
  rm(mask)
  gc()

  print("Created time series to matrix")


  m <- parallel::mclapply(voxNames,
                          FUN = listFormula, formula, mc.cores = ncores)

  print("Created formula list")

  timeIn <- proc.time()
  imageMat <- cbind(imageMat, subjData)

  print("Running test model")
  model <- lmerTest::lmer(m[[1]], data=imageMat, ...)

  print("Running parallel models")
  model <- parallel::mclapply(m,
                              FUN = function(x, data,  ...) {
                                foo <- base::do.call(lmerTest::lmer, list(formula = x, data=data, ...))
                                foo <- base::do.call(lmerTest::lmer, list(formula = x, data=data, ...))
                                return(list(summary(foo)$coefficients, summary(foo)$residuals))
                              }, data=imageMat, ..., mc.preschedule = mc.preschedule, mc.cores = ncores)

  timeOut <- proc.time() - timeIn
  print(timeOut[3])
  print("Parallel Models Ran")

  return(model)

}
