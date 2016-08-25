#' Timeseries to Mean Cluster 
#'
#' This function is able to output the mean voxel intensity over a cluster. 
#' Each column represents a cluster and the rows represent the t-dimention.
#' @param image Input image of type 'nifti'
#' @param mask Input mask of type 'nifti'. Must have different clusters labeled as integers.
#' 
#' @export
#' @examples
#' 
#' 
#' image <- oro.nifti::nifti(img = array(1:320, dim =c(4,4,4,5)))
#' mask <- oro.nifti::nifti(img = array(0:15, dim = c(4,4,4,1)))
#' matrix <- ts2meanCluster(image, mask)


ts2meanCluster <- function(image, mask) {
  
  if (missing(image)) { stop("image is missing")}
  if (missing(mask)) { stop("mask is missing")}
  
  label <- sort(as.numeric(unique(matrix(mask@.Data))))
  label <- label[label != 0]
  
  if (!all.equal(round(label), round(label))) {
    stop("Cluster Mask must have integer labels")
  }
  
  if (!all.equal(dim(image@.Data)[1:3], dim(mask@.Data)[1:3])) {
    stop("Image and Mask have different dimentions")
  }
  
  matrix <- matrix(NA, nrow = dim(image@.Data)[length(dim(image@.Data))], ncol=max(label))
  
  for (i in label) {
    temp <- matrix(image@.Data)[mask@.Data == i]
    dim(temp) <- c(sum(mask@.Data == i), dim(image@.Data)[length(dim(image@.Data))])
    temp <- t(temp)
    matrix[,i] <- rowMeans(temp)
  }
  
  matrix <- as.data.frame(matrix)
  names <- base::lapply(1:dim(matrix)[2], function(x) { return(paste0("cluster",x))})
  names(matrix) <- names
  gc()
  
  keep.NonNA <- base::rep(TRUE, dim(matrix)[2])
  for (i in 1:dim(matrix)[2]) {
    if (all(is.na(matrix[,i]))) {
      keep.NonNA[i] <- FALSE
    }
  }
  
  matrix <- matrix[,keep.NonNA]
  
  return(matrix)
}