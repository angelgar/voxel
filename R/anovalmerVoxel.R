#' Computes voxelwise analysis of variance (ANOVA) tables for a Linear Mixed Effects Model.
#'
#' This function computes analysis of variance tables for the fitted models after running a Linear Mixed Effect Model using the lmerTest() function and the anova function in that package.
#' The analysis will run in all voxels in the mask and will return the analysis of variance table for each voxel.
#' Please check the lmerTest documentation for further information about specific arguments used in lmerTest::anova(). Multi-model calls are disabled.
#'
#'
#'
#' @param image Input image of type 'nifti' or vector of path(s) to images. If multiple paths, the script will call mergeNifti() and merge across time.
#' @param mask Input mask of type 'nifti' or path to mask. Must be a binary mask
#' @param fourdOut To be passed to mergeNifti, This is the path and file name without the suffix to save the fourd file. Default (NULL) means script won't write out 4D image.
#' @param formula Must be a formula passed to lmer()
#' @param subjData Dataframe containing all the covariates used for the analysis
#' @param ddf Which approximation of DDF to be used. To be passed to anova.merModLmerTest. Defaults to "Satterthwaite"
#' @param type Type of hypothesis to be test (defined from SAS theory). Defaults to 3. To be passed to anova.merModLmerTest
#' @param mc.preschedule Argument to be passed to mclapply, whether or not to preschedule the jobs. More info in parallel::mclapply
#' @param ncores Number of cores to use
#' @param ... Additional arguments passed to lmer()
#'
#' @return Returns list of models fitted to each voxel over the masked images passed to function.
#' @export
#'
#'
#'
#' @examples
#'
#'
#' image <- oro.nifti::nifti(img = array(1:1600, dim =c(4,4,4,25)))
#' mask <- oro.nifti::nifti(img = array(c(rep(0,15), rep(1,1)), dim = c(4,4,4,1)))
#' set.seed(1)
#' covs <- data.frame(x = runif(25), y = runif(25), id = rep(1:5,5))
#' fm1 <- "~ x + y + (1|id)"
#' models <- anovalmerVoxel(image, mask, formula = fm1, subjData = covs, ncores = 1, REML=TRUE)
#'


anovalmerVoxel <- function(image, mask , fourdOut = NULL, formula, subjData, ddf="Satterthwaite", type=3, mc.preschedule = TRUE, ncores = 1, ...) {

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

  voxNames <- as.character(names(imageMat))

  rm(image)
  rm(mask)
  gc()

  print("Created time series to matrix")

  m <- parallel::mclapply(voxNames, FUN = listFormula, formula, mc.cores = ncores)

  rm(formula)
  gc()

  imageMat <- cbind(imageMat, subjData)

  print("Created formula list")

  timeIn <- proc.time()
  print("Running test ANOVA")

  foo <- base::do.call(lmerTest::lmer, list(formula = m[[1]], data=imageMat, ...))
  ## foo <- as.call(list(quote(lmerTest::lmer), formula=m[[1]], data=imageMat, ...))
  model <- anova(foo, ddf=ddf, type=type)

  print("Running parallel ANOVAs")
  model <- parallel::mclapply(m,
                              FUN = function(x, data, ddf, type,  ...) {
                                foo <- base::do.call(lmerTest::lmer, list(formula = x, data=data, ...))
                                return(lmerTest_anova(foo, ddf=ddf, type=type))
                              }, data=imageMat,ddf=ddf,type=type, ..., mc.preschedule = mc.preschedule , mc.cores = ncores)

  timeOut <- proc.time() - timeIn
  print(timeOut[3])
  print("Parallel ANOVAs Ran")

  return(model)

}

