#' Generate FSL Randomise call for a GAM Model
#'  
#'
#' This function is able to generate all the necessary files to run randomise with a GAM Model
#' This script will write out all design and contrast files 
#' This function will run a f-test to compare a full and reduced model (a model with and without spline)
#' 
#' 
#' @param image Input path of 'nifti' image or vector of path(s) to images. If multiple paths, the script will all mergeNiftis() and merge across time.
#' @param maskPath to mask. Must be a binary mask
#' @param formulaFull Must be the formula of the full model (i.e. "~s(age,k=5)+sex+mprage_antsCT_vol_TBV")
#' @param formulaRed Must be the formula of the reduced model (i.e. "~sex+mprage_antsCT_vol_TBV")
#' @param subjData Dataframe containing all the covariates used for the analysis
#' @param outDir output directory for randomise 
#' @param nsim Number of simulations
#' @param thresh significance threshold
#' @param run FALSE will only print randomise command but won't it
#' 
#' @return Return randomise command 
#' @export
#' 
#' @examples
#' \dontrun{
#' 
#' subjData = mgcv::gamSim(1,n=400,dist="normal",scale=2)
#' OutDirRoot="Output Directory"
#' maskName="Path to mask"
#' imagePath="Path to output"
#' covsFormula="~s(age,k=5)+sex+mprage_antsCT_vol_TBV"
#' redFormula="~sex+mprage_antsCT_vol_TBV"
#' 
#' gamRandomise(image = imagePath, maskPath = maskName, formulaFull = covsFormula, 
#'            formulaRed = redFormula, subjData = subjData, outDir = OutDirRoot)
#' 
#' }
#' 



gamRandomise <- function(image, maskPath = NULL, formulaFull, formulaRed, 
                     subjData, outDir, 
                     nsim = 500, thresh = 0.01, run = FALSE) {
  
  if (missing(image)) { stop("image is missing")}
  if (missing(formulaFull)) { stop("formulaFull is missing")}
  if (missing(formulaRed)) { stop("formulaRed is missing")}
  if (missing(subjData)) { stop("subjData is missing")}
  if (missing(outDir)) { stop("outDir is missing")}
  
  if (class(formulaFull) != "character") { stop("formula class must be character")}
  if (class(formulaRed) != "character") { stop("formula class must be character")}
  
  if (class(image) == "character" & length(image) == 1) {
    mergednifti <- image
  } else if (class(image) == "character" & length(image) > 1) {
    mergednifti = file.path(outDir, 'fourd')
    image <- mergeNiftis(inputPaths = image, direction = "t", outfile = mergednifti)
    mergednifti <- file.path(outDir, 'fourd.nii.gz')
  }
  
  rm(image)
  
  subjData$dummy <- stats::rnorm(dim(subjData)[1])
  # model matrices
  X = mgcv::model.matrix.gam(mgcv::gam(stats::update.formula(formulaFull, "dummy ~ .") , data=subjData))
  Xred = mgcv::model.matrix.gam(mgcv::gam(stats::update.formula(formulaRed, "dummy ~ .") , data=subjData))
  
  ## DESIGN AND CONTRASTS ##
  # design file
  n = nrow(X)
  p = ncol(X)
  p2 = p - ncol(Xred)
  matfile = file.path(outDir, 'design.mat')
  cat('/NumWaves\t', ncol(X), '\n/NumPoints\t', nrow(X), '\n/PPheights\t', paste(apply(X, 2, function(x) abs(diff(range(x))) ), collapse='\t'), '\n\n/Matrix\n', sep='', file=matfile)
  utils::write.table(X, append=TRUE, file=matfile, row.names=FALSE, col.names=FALSE)
  
  # contrast file
  confile1 = file.path(outDir, 'design.con') # for f-test
  cons = matrix(0, nrow=p2, ncol=ncol(X))
  cons[ cbind(1:(p2), which(! colnames(X) %in% colnames(Xred) ) ) ] = 1
  cat('/ContrastName1\t temp\n/ContrastName2\t\n/NumWaves\t', ncol(X), '\n/NumPoints\t', nrow(cons), '\n/PPheights\t', paste(rep(1,ncol(cons)), collapse='\t'), '\n/RequiredEffect\t1\t1\n\n/Matrix\n', sep='', file=confile1)
  utils::write.table(cons, append=TRUE, file=confile1, row.names=FALSE, col.names=FALSE)
  
  # fts file
  ftsfile = file.path(outDir, 'design.fts')
  fts = matrix(1, nrow=1, ncol=nrow(cons)) # ftest of all contrasts
  cat('/NumWaves\t', nrow(cons), '\n/NumContrasts\t', 1, '\n\n/Matrix\n', sep='', file=ftsfile)
  utils::write.table(fts, append=TRUE, file=ftsfile, row.names=FALSE, col.names=FALSE)
  
  # t distribution is two tailed, F is one tailed. -x outputs voxelwise statistics -N outputs null distribution text files
  # F-test
  
  ##Change mergenifti
  if(!is.null(maskPath)){
    fcmd = paste('randomise -i', mergednifti, '-m', maskPath, '-o', file.path(outDir, 'randomise'), '-d', matfile, '-t', confile1, '-f', ftsfile, '--fonly -F', stats::qf( (1-thresh),df1=p2, df2=(n-p) ), '-x -N -n', nsim, '--uncorrp' )
  } else {
    fcmd = paste('randomise -i', mergednifti, '-o', file.path(outDir, 'randomise'), '-d', matfile, '-t', confile1, '-f', ftsfile, '--fonly -F', stats::qf( (1-thresh),df1=p2, df2=(n-p) ), '-x -N -n', nsim, '--uncorrp' )
  }
  
  if(run){
    system(fcmd)
  }
  
  print(fcmd)
  
}