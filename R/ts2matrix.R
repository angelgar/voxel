#' Timeseries to Matrix 
#'
#' This function is able to mask a 4-Dimentional image and create a matrix from it. 
#' Each column represents the same voxel in the xyz array while the rows represent the t-dimention.
#' @param image Input image of type 'nifti'
#' @param mask Input mask of type 'nifti'. Must be a binary mask
#' @keywords internal
#' @export
#' @examples
#' 
#' 
#' image <- oro.nifti::nifti(img = array(1:64, dim =c(4,4,4,5)))
#' mask <- oro.nifti::nifti(img = array(0:1, dim = c(4,4,4)))
#' matrix <- ts2matrix(image, mask)

ts2matrix <- function(image, mask) {
  
  if (missing(image)) { stop("image is missing")}
  if (missing(mask)) { stop("mask is missing")}
  
  label <- sort(as.numeric(unique(matrix(mask@.Data))))
  
  if (length(label) == 2 && label[1] == 0 && label[2] == 1) {
    if (length(dim(image@.Data)) == 3 | dim(image@.Data)[4] == 1) {
      vector <- image@.Data[mask@.Data == 1]
      gc()
      return(vector)
      
    } else {
      temp <- matrix(image@.Data)[mask@.Data == 1]
      dim(temp) <- c(sum(mask@.Data), dim(image@.Data)[length(dim(image@.Data))])
      temp <- t(temp)
      
      temp <- as.data.frame(temp)
      names <- base::lapply(1:dim(temp)[2], function(x) { return(paste0("voxel",x))})
      names(temp) <- names
      gc()
      return(temp)
    }
    
  } else {
    gc()
    stop("Mask Image is not Binary")
  }
}
