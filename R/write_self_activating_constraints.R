#'\code{write_self_activating_constraints}
#'
#'Loop constraints
#'

write_self_activating_constraints <- function(pknList = pknList, 
                                              binaries = binaries, 
                                              dataMatrix = dataMatrix, 
                                              M = M){
  
  cc1 <- c()
  cc2 <- c()
  cc3 <- c()
  
  sif <- createSIF(pknList = pknList)
  
  species <- unique(c(sif[, 1], sif[, 3]))
  
  distVar <- c()
  for(ii in 1:nrow(dataMatrix[[1]])){
    distVar <- c(distVar, paste0("dist{", species, "}_", ii))
  }
  
  speciesVar <- binaries[[1]][grepl(pattern = "species", x = binaries[[3]])]
  speciesExp <- binaries[[3]][grepl(pattern = "species", x = binaries[[3]])]
  
  for(i in 1:length(distVar)){
    
    meas1 <- strsplit(x = distVar[i], split = "{", fixed = TRUE)[[1]][2]
    meas <- strsplit(x = meas1, split = "}", fixed = TRUE)[[1]][1]
    exp <- 
      strsplit(x = distVar[i], 
               split = "_", 
               fixed = TRUE)[[1]][length(strsplit(x = distVar[i], s
                                                  plit = "_", 
                                                  fixed = TRUE)[[1]])]
    
    cc1 <- c(cc1, 
             paste0(binaries[[1]][which(binaries[[3]]==
                                          paste0("species ", meas, 
                                                 " in experiment ", exp))], 
                    " - ", distVar[i], " <= 0"))
    
  }
  
  ##
  reacVar <- binaries[[1]][grepl(pattern = "interaction", x = binaries[[3]])]
  reacExp <- binaries[[3]][grepl(pattern = "interaction", x = binaries[[3]])]
  cnt <- 1
  for(ii in 1:nrow(dataMatrix[[1]])){
    
    for(i in 1:nrow(sif)){
      
      ss <- sif[i, 1]
      tt <- sif[i, 3]
      
      cc2 <- c(cc2, paste0(paste0("dist{", tt, "}_", ii), " - ", 
                           paste0("dist{", ss, "}_", ii), " - ", M, " ", 
                           reacVar[cnt], " >= ", 1-M))
      
      cnt <- cnt + 1
      
    }
    
  }
  
  ##
  cc3 <- paste0(distVar, " <= ", M)
  
  return(c(cc1, cc2, cc3))
  
}