#' Computes analysis of variance tables for a generalized Additive Model on all voxels of a NIfTI image within a mask.
#'
#' This function computes analysis of variance tables for the fitted models after running a Generalized Additive Model (from mgcv::gam).
#' The analysis will run in all voxels in the mask and will return the analysis of variance table for each voxel.
#' Please check anova::gam to learn more about the function of the arguments that are passed to anova.gam. Multi-model calls are disabled.
#'
#' @param image Input image of type 'nifti' or vector of path(s) to images. If multiple paths, the script will call mergeNifti() and merge across time.
#' @param mask Input mask of type 'nifti' or path to mask. Must be a binary mask
#' @param fourdOut To be passed to mergeNifti, This is the path and file name without the suffix to save the fourd file. Default (NULL) means script won't write out 4D image.
#' @param formula Must be a formula passed to gam()
#' @param subjData Dataframe containing all the covariates used for the analysis
#' @param dispersion To be passed to mgcv::anova.gam, Defaults to NULL. Dispersion Parameter, not normally used.
#' @param freq To be passed to mgcv::anova.gam, Defaults to FALSE. Frequentist or Bayesian approximations for p-values
#' @param p.type To be passed to mgcv::anova.gam, Defaults to 0. Exact test statistics o use for smooth terms.
#' @param mc.preschedule Argument to be passed to mclapply, whether or not to preschedule the jobs. More info in parallel::mclapply
#' @param ncores Number of cores to use
#' @param ... Additional arguments passed to gam()
#'
#' @return List of models fitted to each voxel over the masked images passed to function.
#'
#' @export
#' @examples
#' image <- oro.nifti::nifti(img = array(1:1600, dim =c(4,4,4,25)))
#' mask <- oro.nifti::nifti(img = array(0:1, dim = c(4,4,4,1)))
#' set.seed(1)
#' covs <- data.frame(x = runif(25), y=runif(25))
#' fm1 <- "~ s(x) + y"
#' models <- anovagamVoxel(image=image, mask=mask,
#'               formula=fm1, subjData=covs, ncores = 1)
#' @importFrom mgcv anova.gam gam
anovagamVoxel <- function(image, mask , fourdOut = NULL, formula, subjData, dispersion = NULL,
                          freq = FALSE, p.type=0, mc.preschedule = TRUE, ncores = 1, ...) {

  if (missing(image)) { stop("image is missing")}
  if (missing(mask)) { stop("mask is missing")}
  if (missing(formula)) { stop("formula is missing")}
  if (missing(subjData)) { stop("subjData is missing")}

  if (class(formula) != "character") { stop("formula class must be character")}

  if (class(image) == "character" & length(image) == 1) {
    image <- oro.nifti::readNIfTI(fname=image)
  } else if (class(image) == "character" & length(image) > 1) {
    image <- mergeNiftis(inputPaths = image, direction = "t", outfile <- fourdOut)
  }

  if (class(mask) == "character" & length(mask) == 1) {
    mask <- oro.nifti::readNIfTI(fname=mask)
  }

  imageMat <- ts2matrix(image, mask)

  rm(image)
  rm(mask)
  gc()

  print("Created time series to matrix")

  voxNames <- names(imageMat)

  m <- parallel::mclapply(voxNames,
                          FUN = listFormula, formula, mc.cores = ncores)

  imageMat <- cbind(imageMat, subjData)

  print("Created formula list")

  timeIn <- proc.time()
  print("Running test ANOVA")
  model <- mgcv::gam(m[[1]], data=imageMat, ...)
  model <- mgcv::anova.gam(model, dispersion = dispersion,
                           freq = freq, p.type = p.type)

  print("Running parallel ANOVAs")
  model <- parallel::mclapply(m,
                              FUN = function(x, data, dispersion, freq, p.type, ...) {
                                foo <- mgcv::gam(x, data=data, ...)
                                foo <- mgcv::anova.gam(foo, dispersion = dispersion,
                                                       freq = freq, p.type = p.type)
                                return(list(anovaPTable = foo$pTerms.table, anovaSTable = foo$s.table))
                              }, data=imageMat, dispersion, freq, p.type, ..., mc.preschedule = mc.preschedule, mc.cores = ncores)


  timeOut <- proc.time() - timeIn
  print(timeOut[3])
  print("Parallel ANOVAs Ran")

  return(model)

}
