#' Merge NIfTI Images across specified direction
#'
#' This function merges nifti images together in a specified direction.
#' 
#' @param inputPaths This is a vector of input filenames (character)
#' @param direction This is the direction you want to merge your image over, x, y, z, or t 
#' @param outfile This is the path and file name to save the Nifti file without the suffix, passed to writeNIfTI
#' @param ncores Number of cores to be used for this operation
#' @param ... Additional arguments passed to readNIfTI
#' @return Returns a merged NIfTI image
#' @keywords internal
#' @export



mergeNiftis <- function(inputPaths, direction = c("x","y","z","t"), outfile = NULL, ncores = 1,...) {
  
  if (missing(inputPaths)) { stop("inputPaths is missing")}
  
  if (class(inputPaths) != "character") { stop("inputPaths is not a character vector of paths)")}
  if (length(inputPaths) < 2) { stop("Input Paths has less than two paths")}
  
  images <- parallel::mclapply(inputPaths, 
                               FUN = function(x) {
                                 gc()
                                 foo <- oro.nifti::readNIfTI(fname=x,...)
                                 return(list(foo@.Data, foo@datatype))
                               },..., mc.cores = ncores, mc.preschedule = F)
  
  datatype <- unique(unlist((parallel::mclapply(images, FUN = function(x) {
                                                  gc()
                                                  return(x[[2]])
                                       }, mc.cores=ncores, mc.preschedule = F))))
  
  if (length(unique(datatype)) != 1) {
    stop("Images have different datatypes, cannot merge")
  }
  
  dim <- as.vector(parallel::mclapply(images, FUN = function(x) {
                                                dim(x[[1]])
                                           }, mc.cores=ncores, mc.preschedule = T))
  
  if (length(unique(dim)[[1]]) != length(dim(images[[1]][[1]]))) {
    stop("Images have different number of dimentions, cannot merge")
  }
  
  
  if (any(unique(dim)[[1]] != dim(images[[1]][[1]]))) {
    stop("Images have different dimentions, cannot merge")
  }
  
  if (length(unique(dim)[[1]]) == 4) {
    
    if (direction == "t") {
      
      timeIndex <- length(images) * dim(images[[1]][[1]])[4]
      mergedArray <- array(data=NA, dim=c(dim(images[[1]][[1]])[1:3], timeIndex))
      ImageL <- dim(images[[1]][[1]])[4]
      
      for (i in 1:length(images)) {
        mergedArray[,,,(1 + (i - 1)*ImageL ):((i)*ImageL)] <- images[[i]][[1]]
      }
      
      rm(images)
      
      mergedNifti <- oro.nifti::nifti(img = mergedArray, datatype = datatype)
      rm(mergedArray)
      gc()
    }
    
    
    if (direction == "x") {
      
      timeIndex <- length(images) * dim(images[[1]][[1]])[1]
      mergedArray <- array(data=NA, dim=c( timeIndex, dim(images[[1]][[1]])[2:4]))
      ImageL <- dim(images[[1]][[1]])[1]
      
      for (i in 1:length(images)) {
        mergedArray[(1 + (i - 1)*ImageL ):((i)*ImageL),,,] <- images[[i]][[1]]
      }
      
      rm(images)
      
      mergedNifti <- oro.nifti::nifti(img = mergedArray, datatype = datatype)
      rm(mergedArray)
      gc()
    }
    
    
    if (direction == "y") {
      
      timeIndex <- length(images) * dim(images[[1]][[1]])[2]
      mergedArray <- array(data=NA, dim=c(dim(images[[1]][[1]])[1], timeIndex, dim(images[[1]][[1]])[3:4]))
      ImageL <- dim(images[[1]][[1]])[2]
      
      for (i in 1:length(images)) {
        mergedArray[,(1 + (i - 1)*ImageL ):((i)*ImageL),,] <- images[[i]][[1]]
      }
      
      rm(images)
      
      mergedNifti <- oro.nifti::nifti(img = mergedArray, datatype = datatype)
      rm(mergedArray)
      gc()
    }
    
    
    if (direction == "z") {
      
      timeIndex <- length(images) * dim(images[[1]][[1]])[3]
      mergedArray <- array(data=NA, dim=c(dim(images[[1]][[1]])[1:2], timeIndex, dim(images[[1]][[1]])[4]))
      ImageL <- dim(images[[1]][[1]])[3]
      
      for (i in 1:length(images)) {
        mergedArray[,,(1 + (i - 1)*ImageL ):((i)*ImageL),] <- images[[i]][[1]]
      }
      
      rm(images)
      
      mergedNifti <- oro.nifti::nifti(img = mergedArray, datatype = datatype)
      rm(mergedArray)
      gc()
    }
  }
   
  
  if (length(unique(dim)[[1]]) == 3) {
    
    if (direction == "t") {
      
      timeIndex <- length(images) 
      mergedArray <- array(data=NA, dim=c(dim(images[[1]][[1]])[1:3], timeIndex))
      ImageL <- 1
      
      for (i in 1:length(images)) {
        mergedArray[,,,(1 + (i - 1)*ImageL ):((i)*ImageL)] <- images[[i]][[1]]
      }
      
      rm(images)
      
      mergedNifti <- oro.nifti::nifti(img = mergedArray, datatype = datatype)
      rm(mergedArray)
      gc()
    }
    
    
    if (direction == "x") {
      
      timeIndex <- length(images) * dim(images[[1]][[1]])[1]
      mergedArray <- array(data=NA, dim=c( timeIndex, dim(images[[1]][[1]])[2:3]))
      ImageL <- dim(images[[1]][[1]])[1]
      
      for (i in 1:length(images)) {
        mergedArray[(1 + (i - 1)*ImageL ):((i)*ImageL),,] <- images[[i]][[1]]
      }
      
      rm(images)
      
      mergedNifti <- oro.nifti::nifti(img = mergedArray, datatype = datatype)
      rm(mergedArray)
      gc()
    }
    
    
    if (direction == "y") {
      
      timeIndex <- length(images) * dim(images[[1]][[1]])[2]
      mergedArray <- array(data=NA, dim=c(dim(images[[1]][[1]])[1], timeIndex, dim(images[[1]][[1]])[3]))
      ImageL <- dim(images[[1]][[1]])[2]
      
      for (i in 1:length(images)) {
        mergedArray[,(1 + (i - 1)*ImageL ):((i)*ImageL),] <- images[[i]][[1]]
      }
      
      rm(images)
      
      mergedNifti <- oro.nifti::nifti(img = mergedArray, datatype = datatype)
      rm(mergedArray)
      gc()
    }
    
    
    if (direction == "z") {
      
      timeIndex <- length(images) * dim(images[[1]][[1]])[3]
      mergedArray <- array(data=NA, dim=c(dim(images[[1]][[1]])[1:2], timeIndex))
      ImageL <- dim(images[[1]][[1]])[3]
      
      for (i in 1:length(images)) {
        mergedArray[,,(1 + (i - 1)*ImageL ):((i)*ImageL)] <- images[[i]][[1]]
      }
      
      rm(images)
      
      mergedNifti <- oro.nifti::nifti(img = mergedArray, datatype = datatype)
      rm(mergedArray)
      gc()
    }
  } 
  
  
  if (length(unique(dim)[[1]]) == 2) {
    
    if (direction == "t") {
      stop("Cannot merge a 2D image in the time; chose one of -xyz")
    }
    
    
    if (direction == "x") {
      
      timeIndex <- length(images) * dim(images[[1]][[1]])[1]
      mergedArray <- array(data=NA, dim=c( timeIndex, dim(images[[1]][[1]])[2]))
      ImageL <- dim(images[[1]][[1]])[1]
      
      for (i in 1:length(images)) {
        mergedArray[(1 + (i - 1)*ImageL ):((i)*ImageL),] <- images[[i]][[1]]
      }
      
      rm(images)
      
      mergedNifti <- oro.nifti::nifti(img = mergedArray, datatype = datatype)
      rm(mergedArray)
      gc()
    }
    
    
    if (direction == "y") {
      
      timeIndex <- length(images) * dim(images[[1]][[1]])[2]
      mergedArray <- array(data=NA, dim=c(dim(images[[1]][[1]])[1], timeIndex))
      ImageL <- dim(images[[1]][[1]])[2]
      
      for (i in 1:length(images)) {
        mergedArray[,(1 + (i - 1)*ImageL ):((i)*ImageL)] <- images[[i]][[1]]
      }
      
      rm(images)
      
      mergedNifti <- oro.nifti::nifti(img = mergedArray, datatype = datatype)
      rm(mergedArray)
      gc()
    }
    
    
    if (direction == "z") {
      
      timeIndex <- length(images)
      mergedArray <- array(data=NA, dim=c(dim(images[[1]][[1]])[1:2], timeIndex))
      ImageL <- 1
      
      for (i in 1:length(images)) {
        mergedArray[,,(1 + (i - 1)*ImageL ):((i)*ImageL)] <- images[[i]][[1]]
      }
      
      rm(images)
      
      mergedNifti <- oro.nifti::nifti(img = mergedArray, datatype = datatype)
      rm(mergedArray)
      gc()
    }
  } 
  
  if (!base::is.null(outfile)) {
    oro.nifti::writeNIfTI(mergedNifti, filename=outfile)
  }
  
  
  return(mergedNifti)
  
}