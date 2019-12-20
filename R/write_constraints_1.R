#'\code{write_constraints_1}
#'
#'Wting of constraints - 1
#'

write_constraints_1 <- function(dataMatrix = dataMatrix, 
                                binaries = binaries, pknList = pknList){
  
  sif <- createSIF(pknList)
  
  constraints1 <- c()
  
  for(ii in 1:nrow(sif)){
    
    ss <- sif[ii, 1]
    tt <- sif[ii, 3]
    
    for(jj in 1:nrow(dataMatrix$dataMatrix)){
      
      aa <- 
        binaries[[1]][which(binaries[[3]]==
                              paste0("species ", ss, " in experiment ", jj))]
      bb <- 
        binaries[[1]][which(binaries[[3]]==
                              paste0("species ", tt, " in experiment ", jj))]
      cc <- 
        binaries[[1]][which(binaries[[3]]==
                              paste0("interaction ", ss, "=", tt, 
                                     " in experiment ", jj))]
      
      c1 <- paste0(aa," + ", bb, " - 2", cc, " >= 0")
      
      constraints1 <- c(constraints1, c1)
      
      if(sif[ii, 1]%in%dataMatrix$species[dataMatrix$tgID[[jj]]]){
        
        constraints1 = 
          c(constraints1, paste0(aa, " + ", bb, " - ", cc, " >= 0"))
        constraints1 = 
          c(constraints1, paste0(aa, " + ", bb, " - ", cc, " <= 1"))
        
      }
      
    }
    
  }
  
  constraints1 <- unique(constraints1)
  
  return(constraints1)
  
}
