#'\code{write_boundaries}
#'
#'Set the boundaries for each variables
#'

write_boundaries <- function(binaries = binaries, 
                             pknList = pknList, M = M, dataMatrix = dataMatrix){
  
  bounds <- c()
  for(i in 1:length(binaries[[1]])){
    #bounds <- c(bounds, paste("0 <= ", binaries[[1]][i], " <= 1", sep = ""))
    bounds <- c(bounds, paste("\t", binaries[[1]][i], " >= 0\t \t", sep = ""))
    bounds <- c(bounds, paste("\t", binaries[[1]][i], " <= 1\t \t", sep = ""))
  }
  
  sif <- createSIF(pknList = pknList)
  
  species <- unique(c(sif[, 1], sif[, 3]))
  
  distVar <- c()
  for(ii in 1:nrow(dataMatrix[[1]])){
    distVar <- c(distVar, paste0("dist{", species, "}_", ii))
  }
  
  bounds <- c(bounds, paste0("\t", "0 <= ", distVar, " <= ", M))
  
  return(bounds)
  
}