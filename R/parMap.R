#' Create parametric maps
#'
#' This function create parametric maps according from model parametric tables or analysis of variance tables. 
#' The function will return a p-map, t-map, signed z-map, p-adjusted-map for parametric terms and p-map, z-map, p-adjusted-map for smooth terms. 
#' Additionally the function will return a p-map, F-map, p-toz-map, and p-adjusted-map if the input is ANOVA.
#' You can select which type of p-value correction you want done on the map. The z-maps are signed just like FSL.
#' 
#' @param parameters list of parametric and smooth table coefficents or ANOVA (like the output from vlmParam, vgamParam, anovalmVoxel)
#' @param mask Input mask of type 'nifti' or path to one. Must be a binary mask or a character. Must match the mask passed to one of vlmParam, vgamParam, vgamm4Param, vlmerParam
#' @param method which method of correction for multiple comparisons (default is none)
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
#' models <- vlmParam(image=image, mask=mask, 
#'               formula=fm1, subjData=covs, ncores = 1)
#' Maps <- parMap(models, mask, method="fdr")



parMap <- function(parameters, mask, method="none", outDir = NULL) {
  
  if (missing(parameters)) { stop("parameters is missing")}
  if (missing(mask)) { stop("mask is missing")}
  
  ParameterMaps <- list(NA)
  names(ParameterMaps) <- "temp"
  
  if (class(mask) == "character" & length(mask) == 1) {
    mask <- oro.nifti::readNIfTI(fname=mask)
  }
  
   
  if (class(parameters[[1]])[1] == "anova") {
    if (dim(parameters[[1]])[2] == 5) {
      print("Working with anova.lm object")
      for (j in 1:(dim(parameters[[1]])[1] - 1)) {
        
        pOut<-matrix(NA,nrow=length(parameters),ncol=1)
        FOut<-matrix(NA,nrow=length(parameters),ncol=1)
        zOut<-matrix(NA,nrow=length(parameters),ncol=1)
        
        variable <- rownames(parameters[[1]])[j]
        pvalIndex <- which(colnames(parameters[[1]]) == "Pr(>F)")
        FvalIndex <- which(colnames(parameters[[1]]) == "F value")
        
        for( i in 1:length(parameters)){
          pOut[i,1]<- parameters[[i]][which(rownames(parameters[[i]]) == variable),pvalIndex]
          FOut[i,1]<- parameters[[i]][which(rownames(parameters[[i]]) == variable),FvalIndex]
          zOut[i,1]<- stats::qnorm(parameters[[i]][which(rownames(parameters[[i]]) == variable),pvalIndex] / 2, lower.tail=F)
        }
        
        pOutImage<-mask
        zOutImage<-mask
        FOutImage<-mask
        
        pOutImage@.Data[mask==1@.Data]<-pOut
        zOutImage@.Data[mask==1@.Data]<-zOut
        FOutImage@.Data[mask==1@.Data]<-FOut
        
        pAdjustedOutImage<-mask
        pAdjustedOut <- stats::p.adjust(pOut, method=method)
        pAdjustedOutImage@.Data[mask==1@.Data]<-pAdjustedOut
        
        
        var <- gsub("\\(", "", variable)
        var <- gsub("\\)", "", var)
        var <- gsub("\\*","and",var)
        var <- gsub(":","and",var)
        
        tempNames <- names(ParameterMaps)
        ParameterMaps <- c(ParameterMaps, list(pOutImage), list(FOutImage), list(zOutImage), list(pAdjustedOutImage))
        names(ParameterMaps) <- c(tempNames, paste0(var,"_pAnovaMap"), paste0(var,"_FAnovaMap"), paste0(var, "_zAnovaMap"), paste0(var, "_MultipleComp_pAnovaAdjusted_",method,"_Map"))
        
      }
    } else if (dim(parameters[[1]])[2] == 6) {
      print("Working with anova.merModLmerTest object")
      for (j in 1:(dim(parameters[[1]])[1])) {
        
        pOut<-matrix(NA,nrow=length(parameters),ncol=1)
        FOut<-matrix(NA,nrow=length(parameters),ncol=1)
        zOut<-matrix(NA,nrow=length(parameters),ncol=1)
        
        variable <- rownames(parameters[[1]])[j]
        pvalIndex <- which(colnames(parameters[[1]]) == "Pr(>F)")
        FvalIndex <- which(colnames(parameters[[1]]) == "F.value")
        
        for( i in 1:length(parameters)){
          FOut[i,1]<- parameters[[i]][which(rownames(parameters[[i]]) == variable),FvalIndex]
          pOut[i,1]<- stats::pf(FOut[i,1], parameters[[i]]$NumDF[which(rownames(parameters[[i]]) == variable)],
                         parameters[[i]]$DenDF[which(rownames(parameters[[i]]) == variable)], lower.tail = F)                           
          zOut[i,1]<- stats::qnorm(parameters[[i]][which(rownames(parameters[[i]]) == variable),pvalIndex] / 2, lower.tail=F)
        }
        
        pOutImage<-mask
        zOutImage<-mask
        FOutImage<-mask
        
        pOutImage@.Data[mask==1@.Data]<-pOut
        zOutImage@.Data[mask==1@.Data]<-zOut
        FOutImage@.Data[mask==1@.Data]<-FOut
        
        pAdjustedOutImage<-mask
        pAdjustedOut <- stats::p.adjust(pOut, method=method)
        pAdjustedOutImage@.Data[mask==1@.Data]<-pAdjustedOut
        
        
        var <- gsub("\\(", "", variable)
        var <- gsub("\\)", "", var)
        var <- gsub("\\*","and",var)
        var <- gsub(":","and",var)
        
        tempNames <- names(ParameterMaps)
        ParameterMaps <- c(ParameterMaps, list(pOutImage), list(FOutImage), list(zOutImage), list(pAdjustedOutImage))
        names(ParameterMaps) <- c(tempNames, paste0(var,"_pAnovaMap"), paste0(var,"_FAnovaMap"), paste0(var, "_zAnovaMap"), paste0(var, "_MultipleComp_pAnovaAdjusted_",method,"_Map"))
        
      }
    }
    
  } 
  
  if (!is.null(names(parameters[[1]]))) {
    if (length(names(parameters[[1]])) == 2) {
      if (all(names(parameters[[1]]) == c("anovaPTable","anovaSTable"))) {
        print("Working with anova.gam object")
        
        if (!is.null(parameters[[1]]$anovaPTable)) {
          for (j in 1:(dim(parameters[[1]]$anovaPTable)[1])) {
            
            pOut<-matrix(NA,nrow=length(parameters),ncol=1)
            FOut<-matrix(NA,nrow=length(parameters),ncol=1)
            zOut<-matrix(NA,nrow=length(parameters),ncol=1)
            
            variable <- rownames(parameters[[1]]$anovaPTable)[j]
            pvalIndex <- which(colnames(parameters[[1]]$anovaPTable) == "p-value")
            FvalIndex <- which(colnames(parameters[[1]]$anovaPTable) == "F")
            
            for( i in 1:length(parameters)){
              pOut[i,1]<- parameters[[i]]$anovaPTable[which(rownames(parameters[[i]]$anovaPTable) == variable),pvalIndex]
              FOut[i,1]<- parameters[[i]]$anovaPTable[which(rownames(parameters[[i]]$anovaPTable) == variable),FvalIndex]
              zOut[i,1]<- stats::qnorm(parameters[[i]]$anovaPTable[which(rownames(parameters[[i]]$anovaPTable) == variable),pvalIndex] / 2, lower.tail=F)
            }
            
            pOutImage<-mask
            zOutImage<-mask
            FOutImage<-mask
            
            pOutImage@.Data[mask==1@.Data]<-pOut
            zOutImage@.Data[mask==1@.Data]<-zOut
            FOutImage@.Data[mask==1@.Data]<-FOut
            
            pAdjustedOutImage<-mask
            pAdjustedOut <- stats::p.adjust(pOut, method=method)
            pAdjustedOutImage@.Data[mask==1@.Data]<-pAdjustedOut
            
            
            var <- gsub("\\(", "", variable)
            var <- gsub("\\)", "", var)
            var <- gsub("\\*","and",var)
            var <- gsub(":","and",var)
            
            tempNames <- names(ParameterMaps)
            ParameterMaps <- c(ParameterMaps, list(pOutImage), list(FOutImage), list(zOutImage), list(pAdjustedOutImage))
            names(ParameterMaps) <- c(tempNames, paste0(var,"_pAnovaMap"), paste0(var,"_FAnovaMap"), paste0(var, "_zAnovaMap"), paste0(var, "_MultipleComp_pAnovaAdjusted_",method,"_Map"))
            
          }
        }
        
        
        if (!is.null(parameters[[1]]$anovaSTable)) {
          for (j in 1:(dim(parameters[[1]]$anovaSTable)[1])) {
            
            pOut<-matrix(NA,nrow=length(parameters),ncol=1)
            FOut<-matrix(NA,nrow=length(parameters),ncol=1)
            zOut<-matrix(NA,nrow=length(parameters),ncol=1)
            
            variable <- rownames(parameters[[1]]$anovaSTable)[j]
            pvalIndex <- which(colnames(parameters[[1]]$anovaSTable) == "p-value")
            FvalIndex <- which(colnames(parameters[[1]]$anovaSTable) == "F")
            
            for( i in 1:length(parameters)){
              pOut[i,1]<- parameters[[i]]$anovaSTable[which(rownames(parameters[[i]]$anovaSTable) == variable),pvalIndex]
              FOut[i,1]<- parameters[[i]]$anovaSTable[which(rownames(parameters[[i]]$anovaSTable) == variable),FvalIndex]
              zOut[i,1]<- stats::qnorm(parameters[[i]]$anovaSTable[which(rownames(parameters[[i]]$anovaSTable) == variable),pvalIndex] / 2, lower.tail=F)
            }
            
            pOutImage<-mask
            zOutImage<-mask
            FOutImage<-mask
            
            pOutImage@.Data[mask==1@.Data]<-pOut
            zOutImage@.Data[mask==1@.Data]<-zOut
            FOutImage@.Data[mask==1@.Data]<-FOut
            
            pAdjustedOutImage<-mask
            pAdjustedOut <- stats::p.adjust(pOut, method=method)
            pAdjustedOutImage@.Data[mask==1@.Data]<-pAdjustedOut
            
            
            var <- gsub("\\(", "", variable)
            var <- gsub("\\)", "", var)
            var <- gsub("\\*","and",var)
            var <- gsub(":","and",var)
            
            tempNames <- names(ParameterMaps)
            ParameterMaps <- c(ParameterMaps, list(pOutImage), list(FOutImage), list(zOutImage), list(pAdjustedOutImage))
            names(ParameterMaps) <- c(tempNames, paste0(var,"_pAnovaMap"), paste0(var,"_FAnovaMap"), paste0(var, "_zAnovaMap"), paste0(var, "_MultipleComp_pAnovaAdjusted_",method,"_Map"))
            
          }
        }
      } 
    }
  }
  
  if (length(parameters[[1]]) == 2) {
    if (!all(names(parameters[[1]]) == c("anovaPTable","anovaSTable"))) {
      
      print("Working with output of gam object")
      
      if (!is.null(dim(parameters[[1]][[1]]))) {
        for (j in 1:dim(parameters[[1]][[1]])[1]) {
          
          pOut<-matrix(NA,nrow=length(parameters),ncol=1)
          tOut<-matrix(NA,nrow=length(parameters),ncol=1)
          zOut<-matrix(NA,nrow=length(parameters),ncol=1)
          
          variable <- rownames(parameters[[1]][[1]])[j]
          pvalIndex <- which(colnames(parameters[[1]][[1]]) == "Pr(>|t|)")
          tvalIndex <- which(colnames(parameters[[1]][[1]]) == "t value")
          
          for( i in 1:length(parameters)){
            pOut[i,1]<- parameters[[i]][[1]][which(rownames(parameters[[i]][[1]]) == variable),pvalIndex]
            tOut[i,1]<- parameters[[i]][[1]][which(rownames(parameters[[i]][[1]]) == variable),tvalIndex]
            zOut[i,1]<- sign(tOut[i,1])*stats::qnorm(parameters[[i]][[1]][which(rownames(parameters[[i]][[1]]) == variable),pvalIndex] / 2, lower.tail=F)
          }
          
          pOutImage<-mask
          pAdjustedOutImage<-mask
          zOutImage<-mask
          tOutImage<-mask
          
          pAdjustedOut <- stats::p.adjust(pOut, method=method)
          
          pOutImage@.Data[mask==1@.Data]<-pOut
          pAdjustedOutImage@.Data[mask==1@.Data]<-pAdjustedOut
          zOutImage@.Data[mask==1@.Data]<-zOut
          tOutImage@.Data[mask==1@.Data]<-tOut
          
          var <- gsub("\\(", "", variable)
          var <- gsub("\\)", "", var)
          var <- gsub("\\*","and",var)
          var <- gsub(":","and",var)
          
          tempNames <- names(ParameterMaps)
          ParameterMaps <- c(ParameterMaps, list(pOutImage), list(tOutImage), list(zOutImage), list(pAdjustedOutImage))
          names(ParameterMaps) <- c(tempNames, paste0(var,"_pMap"), paste0(var,"_tMap"), paste0(var, "_zMap"), paste0(var, "_multipleComp_pAdjusted_",method,"_Map"))
        }
      }
      
      if (!is.null(dim(parameters[[1]][[2]]))) {
        for (j in 1:dim(parameters[[1]][[2]])[1]) {
          
          pOut<-matrix(NA,nrow=length(parameters),ncol=1)
          zOut<-matrix(NA,nrow=length(parameters),ncol=1)
          
          variable <- rownames(parameters[[1]][[2]])[j]
          pvalIndex <- which(colnames(parameters[[1]][[2]]) == "p-value")
          
          for(i in 1:length(parameters)){
            pOut[i,1]<- parameters[[i]][[2]][which(rownames(parameters[[i]][[2]]) == variable),pvalIndex]
            zOut[i,1]<- stats::qnorm(parameters[[i]][[2]][which(rownames(parameters[[i]][[2]]) == variable),pvalIndex] / 2, lower.tail=F)
          }
          
          pOutImage<-mask
          pOutImage@.Data[mask@.Data==1]<-pOut
          
          zOutImage<-mask
          zOutImage@.Data[mask@.Data==1]<-zOut
          
          pAdjustedOutImage<-mask
          pAdjustedOut <- stats::p.adjust(pOut, method=method)
          pAdjustedOutImage@.Data[mask==1@.Data]<-pAdjustedOut
          
          var <- gsub("\\(", "", variable)
          var <- gsub("\\)", "", var)
          var <- gsub(",", "", var)
          var <- gsub("=", "", var)
          
          tempNames <- names(ParameterMaps)
          ParameterMaps <- c(ParameterMaps, list(pOutImage), list(zOutImage), list(pAdjustedOutImage))
          names(ParameterMaps) <- c(tempNames, paste0(var,"_pMap"), paste0(var, "_zMap"), paste0(var, "_multipleComp_pAdjusted_",method,"_Map"))
        }
      }
    }
  } 
  
  if (length(colnames(parameters[[1]])) == 4) {
    if (all(colnames(parameters[[1]]) == c("Estimate","Std. Error","t value","Pr(>|t|)"))) {
      print("Working with output from lm model")
      
      for (j in 1:dim(parameters[[1]])[1]) {
        
        pOut<-matrix(NA,nrow=length(parameters),ncol=1)
        tOut<-matrix(NA,nrow=length(parameters),ncol=1)
        zOut<-matrix(NA,nrow=length(parameters),ncol=1)
        
        variable <- rownames(parameters[[1]])[j]
        pvalIndex <- which(colnames(parameters[[1]]) == "Pr(>|t|)")
        tvalIndex <- which(colnames(parameters[[1]]) == "t value")
        
        for( i in 1:length(parameters)){
          pOut[i,1]<- parameters[[i]][which(rownames(parameters[[i]]) == variable),pvalIndex]
          tOut[i,1]<- parameters[[i]][which(rownames(parameters[[i]]) == variable),tvalIndex]
          zOut[i,1]<- sign(tOut[i,1])*stats::qnorm(parameters[[i]][which(rownames(parameters[[i]]) == variable),pvalIndex] / 2, lower.tail=F)
        }
        
        pOutImage<-mask
        zOutImage<-mask
        tOutImage<-mask
        
        pOutImage@.Data[mask==1@.Data]<-pOut
        zOutImage@.Data[mask==1@.Data]<-zOut
        tOutImage@.Data[mask==1@.Data]<-tOut
        
        pAdjustedOutImage<-mask
        pAdjustedOut <- stats::p.adjust(pOut, method=method)
        pAdjustedOutImage@.Data[mask==1@.Data]<-pAdjustedOut
        
        
        var <- gsub("\\(", "", variable)
        var <- gsub("\\)", "", var)
        var <- gsub("\\*","and",var)
        var <- gsub(":","and",var)
        
        tempNames <- names(ParameterMaps)
        ParameterMaps <- c(ParameterMaps, list(pOutImage), list(tOutImage), list(zOutImage), list(pAdjustedOutImage))
        names(ParameterMaps) <- c(tempNames, paste0(var,"_pMap"), paste0(var,"_tMap"), paste0(var, "_zMap"), paste0(var, "_multipleComp_pAdjusted_",method,"_Map"))
      }
    }
  }
  
  if (length(colnames(parameters[[1]])) == 5) {
    
    if (all(colnames(parameters[[1]]) == c("Estimate","Std. Error","df","t value","Pr(>|t|)"))) {
      
      print("Working with output from lmerTest model")
      
      for (j in 1:dim(parameters[[1]])[1]) {
        
        pOut<-matrix(NA,nrow=length(parameters),ncol=1)
        tOut<-matrix(NA,nrow=length(parameters),ncol=1)
        zOut<-matrix(NA,nrow=length(parameters),ncol=1)
        
        variable <- rownames(parameters[[1]])[j]
        pvalIndex <- which(colnames(parameters[[1]]) == "Pr(>|t|)")
        tvalIndex <- which(colnames(parameters[[1]]) == "t value")
        
        for( i in 1:length(parameters)){
          tOut[i,1]<- parameters[[i]][which(rownames(parameters[[i]]) == variable),tvalIndex]
          pOut[i,1]<- stats::pt(abs(tOut[i,1]), parameters[[i]][which(rownames(parameters[[i]]) == variable), 
                                                         which(colnames(parameters[[1]]) == "df")], lower.tail=F)*2
          zOut[i,1]<- sign(tOut[i,1])*stats::qnorm(parameters[[i]][which(rownames(parameters[[i]]) == variable),pvalIndex] / 2, lower.tail=F)
        }
        
        pOutImage<-mask
        zOutImage<-mask
        tOutImage<-mask
        
        pOutImage@.Data[mask==1@.Data]<-pOut
        zOutImage@.Data[mask==1@.Data]<-zOut
        tOutImage@.Data[mask==1@.Data]<-tOut
        
        pAdjustedOutImage<-mask
        pAdjustedOut <- stats::p.adjust(pOut, method=method)
        pAdjustedOutImage@.Data[mask==1@.Data]<-pAdjustedOut
        
        
        var <- gsub("\\(", "", variable)
        var <- gsub("\\)", "", var)
        var <- gsub("\\*","and",var)
        var <- gsub(":","and",var)
        
        tempNames <- names(ParameterMaps)
        ParameterMaps <- c(ParameterMaps, list(pOutImage), list(tOutImage), list(zOutImage), list(pAdjustedOutImage))
        names(ParameterMaps) <- c(tempNames, paste0(var,"_pMap"), paste0(var,"_tMap"), paste0(var, "_zMap"), paste0(var, "_multipleComp_pAdjusted_",method,"_Map"))
      }
    }
  }
  
  ParameterMaps <- ParameterMaps[-1]
  
  if (length(ParameterMaps) == 0) {
    stop("Input object not recognized")
  }
  
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
  
  return(ParameterMaps)
}
