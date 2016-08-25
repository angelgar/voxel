#' Run a Linear Mixed Effet Model on the mean intensity over a region of interest 
#'
#' This function is able to run a LME using the lmer() function. 
#' All clusters or region of interest must be labeled with integers in the mask passed as an argument.
#' The function relies on lmerTest to create p-values using the Satterthwaite Approximation.
#' 
#' 
#' @param image Input image of type 'nifti' or vector of path(s) to images. If multiple paths, the script will call mergeNifti() and merge across time.
#' @param mask Input mask of type 'nifti' or path to mask. All clusters must be labeled with integers in the mask passed as an argument 
#' @param fourdOut To be passed to mergeNifti, This is the path and file name without the suffix to save the fourd file. Default (NULL) means script won't write out 4D image. 
#' @param formula Must be a formula passed to lmer()
#' @param subjData Dataframe containing all the covariates used for the analysis
#' @param mc.preschedule Argument to be passed to mclapply, whether or not to preschedule the jobs. More info in parallel::mclapply
#' @param ncores Number of cores to use
#' @param ... Additional arguments passed to lmer()
#' @export
#' 
#' @return Returns list of models fitted to the mean voxel intensity a region or interest.
#' 
#' @examples
#' 
#' 
#' image <- oro.nifti::nifti(img = array(1:1600, dim =c(4,4,4,25)))
#' mask <- oro.nifti::nifti(img = array(c(rep(0,14),1,2), dim = c(4,4,4,1)))
#' set.seed(1)
#' covs <- data.frame(x = runif(25), id = rep(1:5,5))
#' fm1 <- "~ x + (1|id)"
#' models <- lmerCluster(image, mask, formula = fm1, subjData = covs, ncores = 1, REML=TRUE)
#' 


lmerCluster <- function(image, mask , fourdOut = NULL, formula, subjData, mc.preschedule = TRUE, ncores = 1, ...) {
  
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
  
  imageMat <- ts2meanCluster(image, mask)
  
  voxNames <- names(imageMat)
  
  rm(image)
  rm(mask)
  gc()
  
  print("Created meanCluster Matrix")
  
  m <- parallel::mclapply(voxNames, 
                          FUN = listFormula, formula, mc.cores = ncores)
  
  imageMat <- cbind(imageMat, subjData) 
  
  print("Created formula list")
  
  timeIn <- proc.time()
  print("Running test model")
  model <- lmerTest::lmer(m[[1]], data=imageMat, ...)
  
  print("Running parallel models")
  model <- parallel::mclapply(m, 
                              FUN = function(x, data,  ...) {
                                base::do.call(lmerTest::lmer, list(formula=x, data=data, ...))
                              }, data=imageMat, ..., mc.preschedule = mc.preschedule , mc.cores = ncores)
  
  timeOut <- proc.time() - timeIn
  print(timeOut[3])
  print("Parallel Models Ran")
  
  return(model)
  
}