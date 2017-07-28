##
buildDataMatrix2 <- function(dataMatrix, sif, targets, experiments){
  
  speciesInNetwork <- unique(c(sif[, 1], sif[, 3]))
  tgID <- which(speciesInNetwork%in%dataMatrix$species[dataMatrix$tgID])
  dnID <- which(speciesInNetwork%in%dataMatrix$species[dataMatrix$dnID])
  dsID <- which(speciesInNetwork%in%dataMatrix$species[dataMatrix$dsID])
  
  dataMatrix2 <- list()
  dataMatrix2$species <- c("")
  for(i in 1:length(tgID)){
    
    dataMatrix2[[1]] <- c(dataMatrix2[[1]], speciesInNetwork[tgID[i]])
    
  }
  for(i in 1:length(dnID)){
    
    dataMatrix2[[1]] <- c(dataMatrix2[[1]], speciesInNetwork[dnID[i]])
    
  }
  for(i in 1:length(dsID)){
    
    dataMatrix2[[1]] <- c(dataMatrix2[[1]], speciesInNetwork[dsID[i]])
    
  }
  dataMatrix2[[1]] <- dataMatrix2[[1]][2:length(dataMatrix2[[1]])]
  
  dataMatrix2$tgID <- 1:length(tgID)
  dataMatrix2$dnID <- (length(tgID)+1):(length(dnID)+length(tgID))
  dataMatrix2$dsID <- (length(tgID)+length(dnID)+1):(length(dsID)+length(dnID)+length(tgID))
  
  return(dataMatrix2)
  
}

##
create_binary_variables_for_x_vector_step_2 <- function(dataMatrix = dataMatrix){
  
  binX <- list()
  
  binX[[1]] <- 1:length(dataMatrix$species)
  
  binX[[2]] <- c("")
  for(i in 2:length(dataMatrix$species)){
    
    binX[[2]] <- c(binX[[2]], paste("x_", binX[[1]][i], sep = ""))
    
  }
  binX[[2]] <- binX[[2]][2:length(binX[[2]])]
  
  binX[[3]] <- c("")
  for(i in 1:length(dataMatrix$species)){
    
    binX[[3]] <- c(binX[[3]], paste("species", dataMatrix$species[i]))
    
  }
  binX[[3]] <- binX[[3]][2:length(binX[[3]])]
  
  return(binX)
  
}

##
create_binary_variables_for_y_vector_step_2 <- function(sif = sif){
  
  y_vector <- c()
  for(i in 1:nrow(sif)){
    y_vector <- c(y_vector, paste(sif[i, 1], "=", sif[i, 3], sep = ""))
  }
  
  variables = c()
  for(i in 1:length(y_vector)){
    variables[i] = paste0("y_",i)
  }
  binary_variables_list = list(1:(length(y_vector)),variables, paste("reaction",y_vector))
  
  return(binary_variables_list)
  
}

##
write_objective_function_step_2 <- function(binaries = binaries, dataMatrix = dataMatrix){
  
  objectiveFunction <- "obj:\t"
  for(i in 1:length(binaries[[3]])){
    
    if(strsplit(binaries[[3]][i], split = " ")[[1]][1] == "species"){
      
      if(!(strsplit(binaries[[3]][i], split = " ")[[1]][2] %in% dataMatrix$species[dataMatrix$dsID])){
        
        objectiveFunction <- paste(objectiveFunction, " + ", binaries[[1]][i], sep = "")
        
      }
      
    }
    
  }
  
  return(objectiveFunction)
  
}

##
write_boundaries_step_2 <- function(binaries = binaries){
  
  bounds <- c()
  for(i in 1:length(binaries[[1]])){
    #bounds <- c(bounds, paste("0 <= ", binaries[[1]][i], " <= 1", sep = ""))
    bounds <- c(bounds, paste("\t", binaries[[1]][i], " >= 0\t \t", sep = ""))
    bounds <- c(bounds, paste("\t", binaries[[1]][i], " <= 1\t \t", sep = ""))
  }
  
  return(bounds)
  
}

##
write_fixed_nodes_constraints <- function(dataMatrix = dataMatrix, binaries = binaries){
  
  fixedConstraints <- c()
  
  for(i in 1:length(binaries[[3]])){
    
    if(strsplit(binaries[[3]][i], split = " ")[[1]][1] == "species"){
      
      if(strsplit(binaries[[3]][i], split = " ")[[1]][2] %in% dataMatrix$species[dataMatrix$tgID]){
        
        fixedConstraints <- c(fixedConstraints, paste(binaries[[1]][i], " = 1", sep = ""))
        
      }
      
      if(strsplit(binaries[[3]][i], split = " ")[[1]][2] %in% dataMatrix$species[dataMatrix$dsID]){
        
        fixedConstraints <- c(fixedConstraints, paste(binaries[[1]][i], " = 1", sep = ""))
        
      }
      
    }
    
  }
  
  return(fixedConstraints)
  
}

##
write_one_edge_constraint <- function(sif = sif, dataMatrix = dataMatrix, binaries = binaries){
  
  reactionPool <- matrix(0, nrow = 1, ncol = 2)
  for(i in 1:length(binaries[[3]])){
    
    if(strsplit(binaries[[3]][i], split = " ")[[1]][1] == "reaction"){
      
      reactionPool <- rbind(reactionPool, c(strsplit(binaries[[3]][i], split = " ")[[1]][2], as.numeric(i)))
      
    }
  
  }
  reactionPool <- reactionPool[2:nrow(reactionPool), ]
  
  constraint <- c()
  
  for(i in 1:length(dataMatrix$dsID)){
    
    kk <- which(sif[, 3]==dataMatrix$species[dataMatrix$dsID[i]])
    currConstraint <- ""
    for(j in 1:length(kk)){
      
      if(j==1){
        
        currConstraint <- paste("xb", as.numeric(reactionPool[which(reactionPool[, 1]==paste(sif[kk[j], 1], "=", sif[kk[j], 3], sep = "")), 2]), sep = "")
        
      }
      else{
        
        currConstraint <- paste(currConstraint, " + xb", as.numeric(reactionPool[which(reactionPool[, 1]==paste(sif[kk[j], 1], "=", sif[kk[j], 3], sep = "")), 2]), sep = "")
        
      }
    }
    
    currConstraint <- paste(currConstraint, " >= 1", sep = "")
    constraint <- c(constraint, currConstraint)

  }
  
  return(constraint)
  
}

##
write_intermediate_node_constraints_out <- function(sif = sif, dataMatrix = dataMatrix, binaries = binaries){
  
  nNames <- dataMatrix$species[dataMatrix$dnID]
  
  reactionPool <- matrix(0, nrow = 1, ncol = 2)
  for(i in 1:length(binaries[[3]])){
    
    if(strsplit(binaries[[3]][i], split = " ")[[1]][1] == "reaction"){
      
      reactionPool <- rbind(reactionPool, c(strsplit(binaries[[3]][i], split = " ")[[1]][2], as.numeric(i)))
      
    }
    
  }
  reactionPool <- reactionPool[2:nrow(reactionPool), ]
  
  nAdjacent <- list()
  for(i in 1:length(nNames)){
    
    temp <- c()
    nAdjacent[[length(nAdjacent)+1]] <- c(temp, which(sif[, 1]==nNames[i]))
    
  }
  
  constraints <- c()
  for(i in 1:length(nAdjacent)){
    
    if((length(nAdjacent[[i]])==0) == FALSE){
      
      currConstraint <- ""
      for(j in 1:length(nAdjacent[[i]])){
        
        if(j==1){
          
          currConstraint <- paste(currConstraint, "xb", as.numeric(reactionPool[which(reactionPool[, 1]==paste(sif[which(sif[, 1]==nNames[i])[j], 1], "=", sif[which(sif[, 1]==nNames[i])[j], 3], sep = "")), 2]), sep = "")
          
        }
        else{
          
          currConstraint <- paste(currConstraint, " + xb", as.numeric(reactionPool[which(reactionPool[, 1]==paste(sif[which(sif[, 1]==nNames[i])[j], 1], "=", sif[which(sif[, 1]==nNames[i])[j], 3], sep = "")), 2]), sep = "")
        }
        
      }
      
      currConstraint <- paste(currConstraint, " - ", binaries[[1]][which(binaries[[3]]==paste("species ", nNames[i], sep = ""))], " >= 0", sep = "")
      
      constraints <- c(constraints, currConstraint)
      
    }
    
  }
  
  return(constraints)
  
}

##
write_intermediate_node_constraints_in <- function(sif = sif, dataMatrix = dataMatrix, binaries = binaries){
  
  nNames <- dataMatrix$species[dataMatrix$dnID]
  
  reactionPool <- matrix(0, nrow = 1, ncol = 2)
  for(i in 1:length(binaries[[3]])){
    
    if(strsplit(binaries[[3]][i], split = " ")[[1]][1] == "reaction"){
      
      reactionPool <- rbind(reactionPool, c(strsplit(binaries[[3]][i], split = " ")[[1]][2], as.numeric(i)))
      
    }
    
  }
  reactionPool <- reactionPool[2:nrow(reactionPool), ]
  
  nAdjacent <- list()
  for(i in 1:length(nNames)){
    
    temp <- c()
    nAdjacent[[length(nAdjacent)+1]] <- c(temp, which(sif[, 3]==nNames[i]))
    
  }
  
  constraints <- c()
  for(i in 1:length(nAdjacent)){
    
    if((length(nAdjacent[[i]])==0) == FALSE){
      
      currConstraint <- ""
      for(j in 1:length(nAdjacent[[i]])){
        
        if(j==1){
          
          currConstraint <- paste(currConstraint, "xb", as.numeric(reactionPool[which(reactionPool[, 1]==paste(sif[which(sif[, 3]==nNames[i])[j], 1], "=", sif[which(sif[, 3]==nNames[i])[j], 3], sep = "")), 2]), sep = "")
          
        }
        else{
          
          currConstraint <- paste(currConstraint, " + xb", as.numeric(reactionPool[which(reactionPool[, 1]==paste(sif[which(sif[, 3]==nNames[i])[j], 1], "=", sif[which(sif[, 3]==nNames[i])[j], 3], sep = "")), 2]), sep = "")
        }
        
      }
      
      currConstraint <- paste(currConstraint, " - ", binaries[[1]][which(binaries[[3]]==paste("species ", nNames[i], sep = ""))], " >= 0", sep = "")
      
      constraints <- c(constraints, currConstraint)
      
    }
    
  }
  
  return(constraints)
  
}

##
write_reaction_constraints <- function(dataMatrix = dataMatrix, sif = sif, binaries = binaries){
  
  constraints <- c()
  
  reactionPool <- matrix(0, nrow = 1, ncol = 2)
  for(i in 1:length(binaries[[3]])){
    
    if(strsplit(binaries[[3]][i], split = " ")[[1]][1] == "reaction"){
      
      reactionPool <- rbind(reactionPool, c(strsplit(binaries[[3]][i], split = " ")[[1]][2], as.numeric(i)))
      
    }
    
  }
  reactionPool <- reactionPool[2:nrow(reactionPool), ]
  
  for(i in 1:nrow(reactionPool)){
    
    ch1 <- binaries[[1]][which(binaries[[3]]==paste("species ", strsplit(reactionPool[i, 1], split = "=")[[1]][1], sep = ""))]
    ch2 <- binaries[[1]][which(binaries[[3]]==paste("species ", strsplit(reactionPool[i, 1], split = "=")[[1]][2], sep = ""))]
    ch3 <- binaries[[1]][which(binaries[[3]]==paste("reaction ", strsplit(reactionPool[i, 1], split = "=")[[1]][1], "=", strsplit(reactionPool[i, 1], split = "=")[[1]][2], sep = ""))]
    
    constraints <- c(constraints, paste(ch1, " + ", ch2, " - 2 ", ch3, " >= 0", sep = ""))
    
    #constraints <- c(constraints, paste(ch1, " + ", ch2, " - ", ch3, " <= 1", sep = ""))
    
  }
  
  return(constraints)
}

##
all_constraints_2 <- function(constraints1 = constraints1, constraints2 = constraints2, constraints3 = constraints3, constraints4 = constraints4, constraints5 = constraints5){
  
  kk <- 1
  allConstraints <- c()
  
  for(i in 1:length(constraints1)){
    
    allConstraints <- c(allConstraints, paste("c", kk, ":\t", constraints1[i], "\t \t", sep = ""))
    kk <- kk + 1
    
  }
  
  for(i in 1:length(constraints2)){
    
    allConstraints <- c(allConstraints, paste("c", kk, ":\t", constraints2[i], "\t \t", sep = ""))
    kk <- kk + 1
    
  }
  
  for(i in 1:length(constraints3)){
    
    allConstraints <- c(allConstraints, paste("c", kk, ":\t", constraints3[i], "\t \t", sep = ""))
    kk <- kk + 1
    
  }
  
  for(i in 1:length(constraints4)){
    
    allConstraints <- c(allConstraints, paste("c", kk, ":\t", constraints4[i], "\t \t", sep = ""))
    kk <- kk + 1
    
  }
  
  for(i in 1:length(constraints5)){
    
    allConstraints <- c(allConstraints, paste("c", kk, ":\t", constraints5[i], "\t \t", sep = ""))
    kk <- kk + 1
    
  }
  
  return(allConstraints)
  
}

##
readOutSIF_step_2<- function(cplexSolutionFileName, binaries = binaries){
  
  reacIDX <- c()
  for(i in 1:length(binaries[[3]])){
    if(strsplit(binaries[[3]][i], split = " ")[[1]][1] == "reaction"){
      reacIDX <- c(reacIDX, i)
    }
  }
  reacVar <- binaries[[1]][reacIDX]
  
  cplexSolutionData <- xmlParse(cplexSolutionFileName)
  cplexSolution <- xmlToList(cplexSolutionData)
  
  cplexSolutionEdges <- list()
  for(i in 1:length(cplexSolution$variables)){
    
    if(cplexSolution$variables[i]$variable[1] %in% reacVar){
      
      cplexSolutionEdges[[length(cplexSolutionEdges)+1]] <- cplexSolution$variables[i]
      
    }
    
  }
  
  reac <- c()
  for(i in 1:length(cplexSolutionEdges)){
    
    if((as.numeric(cplexSolutionEdges[[i]]$variable[3])) >= 0.5){
      
      reacTemp <- binaries[[3]][which(binaries[[1]]==cplexSolutionEdges[[i]]$variable[1])]
      reac <- c(reac, reacTemp)
      #reacTemp <-  binaries[[1]][which(binaries[[1]]==cplexSolutionEdges[[i]]$variable[1])]
      #reac <- c(reac, reacTemp)
      
    }
    
  }
  
  sif <- matrix(, nrow = length(reac), ncol = 3)
  sif[, 2] <- 1
  for(i in 1:length(reac)){
    
    sif[i, 1] <- strsplit(strsplit(reac[i], split = " ")[[1]][2], split = "=")[[1]][1]
    sif[i, 3] <- strsplit(strsplit(reac[i], split = " ")[[1]][2], split = "=")[[1]][2]
    
  }
  
  return(sif)
  
}

##
writeFile <- function(objectiveFunction, constraints, bounds, binaries){
  
  data2 = "testFile2.lp"
  write("enter Problem", data2)
  write("", data2, append = TRUE)
  write("Minimize", data2, append = TRUE)
  
  write(objectiveFunction, data2, append = TRUE)
  
  write("Subject To", data, append = TRUE)
  write(constraints, data, append = TRUE)
  
  write("End", data2, append = TRUE)
  
}

