#' Computes voxelwise analysis of variance (ANOVA) tables for a Generalized Additive Mixed Effects Model. 
#'
#' This function computes analysis of variance tables for the fitted Generalized Additive Mixed Effects (from gamm4::gamm4) models. 
#' The analysis will run in all voxels in the specified mask and will return a list with the ANOVA table at each voxel.
#' Please check the mgcv::anova.gam documentation for further information about specific arguments used in anova.gam. Multi-model calls are disabled.
#' 
#' 
#' @param image Input image of type 'nifti' or vector of path(s) to images. If multiple paths, the script will call mergeNifti() and merge across time.
#' @param mask Input mask of type 'nifti' or path to mask. Must be a binary.
#' @param fourdOut To be passed to mergeNifti, This is the output path to write out the fourd file. Do not include a suffix (i.e. .nii.gz). Default (NULL) means script won't write out 4D image.
#' @param formula Must be a formula passed to gamm4()
#' @param randomFormula Random effects formula passed to gamm4()
#' @param subjData Dataframe containing all the covariates used for the analysis
#' @param dispersion To be passed to mgcv::anova.gam, Defaults to NULL. Dispersion Parameter, not normally used.
#' @param freq To be passed to mgcv::anova.gam, Defaults to FALSE. Frequentist or Bayesian approximations for p-values
#' @param p.type To be passed to mgcv::anova.gam, Defaults to 0. Exact test statistics o use for smooth terms.
#' @param mc.preschedule To be passed to mclapply, whether or not to preschedule the jobs. More info in parallel::mclapply
#' @param ncores Number of cores to use
#' @param ... Additional arguments passed to gamm4()
#' 
#' @return Returns list of models fitted to each voxel over the masked images passed to function.
#' @export
#' 
#' @examples
#' 
#' 
#' image <- oro.nifti::nifti(img = array(1:1600, dim =c(4,4,4,25)))
#' mask <- oro.nifti::nifti(img = array(data = c(rep(0,15), rep(1,1)), 
#'                                      dim = c(4,4,4,1)))
#' set.seed(1)
#' covs <- data.frame(x = runif(25), y=runif(25), id = rep(1:5,5))
#' f1 <- "~ s(x) + y"
#' randomFormula <- "~(1|id)"
#' models <- anovagammVoxel(image, mask, formula = f1, 
#'                               randomFormula = randomFormula, 
#'                               subjData = covs, ncores = 1, REML=TRUE)


anovagammVoxel <- function(image, mask , fourdOut = NULL, formula, randomFormula, 
                                 subjData, dispersion = NULL, freq = FALSE, 
                                 p.type=0, mc.preschedule = TRUE, ncores =1, ...) {
  
  if (missing(image)) { stop("image is missing")}
  if (missing(mask)) { stop("mask is missing")}
  if (missing(formula)) { stop("formula is missing")}
  if (missing(subjData)) { stop("subjData is missing")}
  if (missing(randomFormula)) { stop("randomFormula is missing")}
  
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
  
  print("Created time series to matrix")
  
  rm(image)
  rm(mask)
  gc()
  
  
  voxNames <- names(imageMat)
  
  m <- parallel::mclapply(voxNames, 
                          FUN = listFormula, formula, mc.cores = ncores)
  
  gc()
  
  print("Created formula list")
  
  
  imageMat <- cbind(imageMat, subjData) 
  
  random <- stats::as.formula(randomFormula, env = environment())
  
  timeIn <- proc.time()
  print("Running test ANOVA")
  model <- gamm4::gamm4(m[[1]], data=imageMat, random = random, ...)
  model <- mgcv::anova.gam(model$gam, dispersion = dispersion, 
                           freq = freq, p.type = p.type)
  
  gc()
  print("Running parallel ANOVAs")
  model <- parallel::mclapply(m, 
                              FUN = function(x, data, random, dispersion, freq, p.type,  ...) {
                                foo <- gamm4::gamm4(x, data=data, random=random, ...)
                                foo <- mgcv::anova.gam(foo$gam, dispersion = dispersion, 
                                                       freq = freq, p.type = p.type)
                                return(list(anovaPTable = foo$pTerms.table, anovaSTable = foo$s.table))
                              }, data=imageMat,random=random, dispersion, freq, p.type, ...,mc.preschedule=mc.preschedule, mc.cores = ncores)
  timeOut <- proc.time() - timeIn
  print(timeOut[3])
  
  print("Parallel ANOVAs Ran")
  return(model)
  
}