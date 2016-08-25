#' Run a Generalized Additive Mixed Effects Model on all voxels of a NIfTI image within a mask.  
#'
#' This function is able to run a Generalized Mixed Effects Model (GAMM) using the gamm4() function. 
#' The analysis will run in all voxels within the mask and will return the model fit for each voxel.
#' 
#' 
#' @param image Input image of type 'nifti' or vector of path(s) to images. If multiple paths, the script will call mergeNifti() and merge across time.
#' @param mask Input mask of type 'nifti' or path to mask. Must be a binary mask
#' @param fourdOut To be passed to mergeNifti, This is the path and file name without the suffix to save the fourd file. Default (NULL) means script won't write out 4D image.
#' @param formula Must be a formula passed to gamm4()
#' @param randomFormula Random effects formual passed to gamm4()
#' @param subjData Dataframe containing all the covariates used for the analysis
#' @param mc.preschedule Argument to be passed to mclapply, whether or not to preschedule the jobs. More info in parallel::mclapply
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
#' mask <- oro.nifti::nifti(img = array(c(rep(0,14),1), dim = c(4,4,4,1)))
#' set.seed(1)
#' covs <- data.frame(x = runif(25), id = rep(1:5,5))
#' fm1 <- "~ s(x)"
#' randomFormula <- "~(1|id)"
#' models <- gammVoxel(image = image , mask = mask, formula = fm1, randomFormula = randomFormula, 
#'                                                        subjData = covs, ncores = 1, REML=TRUE)
#' 


gammVoxel <- function(image, mask , fourdOut = NULL, formula, randomFormula, subjData, mc.preschedule = TRUE, ncores =1, ...) {
  
  if (missing(image)) { stop("image is missing")}
  if (missing(mask)) { stop("mask is missing")}
  if (missing(formula)) { stop("formula is missing")}
  if (missing(subjData)) { stop("subjData is missing")}
  if (missing(randomFormula)) { stop("randomFormula is missing")}
  
  if (class(formula) != "character") { stop("formula class must be character")}
  if (class(randomFormula) != "character") { stop("randomFormula class must be character")}
  
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
  print("Running test model")
  model <- gamm4::gamm4(m[[1]], data=imageMat, random = random, ...)
  
  gc()
  print("Running parallel models")
  model <- parallel::mclapply(m, 
                              FUN = function(x, data, random,  ...) {
                                gamm4::gamm4(x, data=data, random=random, ...)
                              }, data=imageMat,random=random, ...,mc.preschedule=mc.preschedule, mc.cores = ncores)
  timeOut <- proc.time() - timeIn
  print(timeOut[3])
  
  print("Parallel Models Ran")
  return(model)
  
}