#' Create parametric maps and residuals
#'
#' This function create parametric maps according from model parametric tables or analysis of variance tables. 
#' The function will return a p-map, t-map, signed z-map, p-adjusted-map for parametric terms and p-map, z-map, p-adjusted-map for smooth terms. 
#' Additionally the function will return a p-map, F-map, p-toz-map, and p-adjusted-map if the input is ANOVA.
#' This function will return a residual map that can be used for cluster correction
#' You can select which type of p-value correction you want done on the map. The z-maps are signed just like FSL.
#' 
#' @param parameters list of parametric and smooth table coefficents or ANOVA (like the output from vlmParam, vgamParam, anovalmVoxel)
#' @param image Input image of type 'nifti' or vector of path(s) to images. If multiple paths, the script will all mergeNifti() and merge across time.
#' @param mask Input mask of type 'nifti' or path to one. Must be a binary mask or a character. Must match the mask passed to one of vlmParam, vgamParam, vgamm4Param, vlmerParam
#' @param method which method of correction for multiple comparisons (default is none)
#' @param ncores Number of cores to use
#' @param mc.preschedule Argument to be passed to mclapply, whether or not to preschedule the jobs. More info in parallel::mclapply
#' @param outDir Path to the folder where to output parametric maps (Default is Null, only change if you want to write maps out)
#' 
#' @return Return parametric maps of the fitted models
#' 
#' @export
#' @examples
#' image <- oro.nifti::nifti(img = array(1:1600, dim =c(4,4,4,25)))
#' mask <- oro.nifti::nifti(img = array(0:1, dim = c(4,4,4,1)))
#' set.seed(1)
#' covs <- data.frame(x = runif(25), y = runif(25))
#' fm1 <- "~ x + y"
#' models <- rlmParam(image=image, mask=mask, 
#'               formula=fm1, subjData=covs, ncores = 1)
#' Maps <- rparMap(models, mask, method="fdr")




rparMap <- function(parameters, image, mask, method, ncores, mc.preschedule, outDir = NULL) {
 
  #Generate tsresiduals
  residualList <- mclapply(parameters, function(x) {
    return(x[[2]])
  }, mc.cores = ncores)
  
  #Generate tsresiduals
  residualMat <- mcmapply(function(x) {
    return(x)
  }, residualList, mc.cores = ncores, SIMPLIFY = TRUE)
  
  rm(residualList)
  gc()
  
  #Save only parameter tables under models
  parameters <- mclapply(parameters, function(x) {
    return(x[[1]])
  }, mc.cores = ncores)
  
  ### Create output
  residualMask <- mask
  residualMask <- residualMask@.Data
  
  #remove image in for memorize optimization purposes
  dataTypeIn <- datatype(image)
  dimPixIn <- pixdim(image)
  rm(image)
  gc()
  
  seq <- 1:dim(residualMat)[1]
  
  #generate 4d residual image
  residuals <- mcmapply(function(x) {
    residualMask[mask@.Data==1] <- residualMat[x,] 
    return(residualMask)
  }, seq, SIMPLIFY = "array", mc.cores = ncores, mc.preschedule= mc.preschedule)
  
  
  #Write it out 
  residualNii <- nifti(residuals, datatype=dataTypeIn, pixdim=dimPixIn)
  rm(residuals)
  gc()
  
  ParameterMaps <- parMap(parameters, mask, method=method)
  ParameterMaps$residuals <- residualNii
  
  rm(parameters)
  gc()
  
  if (!is.null(outDir)) {
    
    dirPath <- base::paste(strsplit(outDir, "/")[[1]][1:(length(strsplit(outDir, "/")[[1]]))], collapse = "/")
    print(base::paste("Checking if", dirPath ,"Exists"))
    if (!dir.exists(dirPath)) {
      print("Directory is missing, creating it now")
      dir.create(dirPath)
    }
    
    for (i in 1:length(names(ParameterMaps))) {
      outPath <- base::paste(outDir, names(ParameterMaps)[i], sep = "/")
      oro.nifti::writeNIfTI(ParameterMaps[[i]], filename = outPath, 
                            gzipped = T)
    }
  }
   
}