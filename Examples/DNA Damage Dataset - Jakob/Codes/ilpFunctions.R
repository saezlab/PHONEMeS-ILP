createSIF <- function(pknList = pknList){
  
  if((class(pknList)=="PKNlist")==FALSE){
    stop("Input object not of class PKNlist")
  }
  
  allInteractions <- pknList@interactions
  SIF <- matrix(, nrow = nrow(allInteractions), ncol = 3)
  
  SIF[, 1] <- allInteractions$K.ID
  
  SIF[, 2] <- rep(1, nrow(allInteractions))
  
  SIF[, 3] <- allInteractions$S.cc
  
  return(SIF)
  
}

##
create_binary_variables_for_y_vector <- function(pknList = pknList){
  
  sif <- createSIF(pknList = pknList)
  
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
create_binary_variables_for_x_vector <- function(dataMatrix = dataMatrix){
  
  allSpecies = c()
  for(i in 1:ncol(dataMatrix[[1]])){
    allSpecies <- c(allSpecies, strsplit(colnames(dataMatrix[[1]])[i], ":")[[1]][2])
  }
  
  variables = c()
  for(i in 1:nrow(dataMatrix[[1]])){
    for(j in 1:length(allSpecies)){
      variables <- c(variables, paste("x_", j, "^", i, sep = ""))
    }
  }
  
  identifiers = c()
  for(i in 1:nrow(dataMatrix[[1]])){
    for(j in 1:length(allSpecies)){
      identifiers <- c(identifiers, paste("species", allSpecies[j], "in experiment", i))
    }
  }
  
  binary_variables_list = list(1:length(variables),variables, identifiers)
  
  return(binary_variables_list)
  
}

##
create_binaries <- function(binaries_x = binaries_x, binaries_y = binaries_y){
  
  numbers <- c()
  bins <- append(binaries_x[[1]], binaries_y[[1]]+length(binaries_x[[1]]))
  for(i in 1:length(bins)){
    numbers <- c(numbers, paste("xb", bins[i], sep = ""))
  }
  variables <- append(binaries_x[[2]], binaries_y[[2]])
  identifiers <- append(binaries_x[[3]], binaries_y[[3]])
  
  binaries = list(numbers, variables, identifiers)
  
  return(binaries)
  
}

##
write_objective_function <- function(dataMatrix = dataMatrix, binaries = binaries){
  
  dM <- dataMatrix[[1]]
  kk <- 1
  objectiveFunction <- "obj:\t"
  for(i in 1:nrow(dM)){
    for(j in 1:ncol(dM)){
      if(as.numeric(dM[i, j]) <= 0){
        if(as.numeric(dM[i, j]) == 0){
          objectiveFunction <- paste(objectiveFunction, "", sep = "")
          kk <- kk+1
        }
        else{
          objectiveFunction <- paste(objectiveFunction, " ", as.numeric(dM[i, j]), " ", "xb", kk, sep = "")
          kk <- kk+1
        }
      }
      else{
        objectiveFunction <- paste(objectiveFunction, " + ", as.numeric(dM[i, j]), " ", "xb", kk, sep = "")
        kk <- kk+1
      }
    }
  }
  
  return(objectiveFunction)
  
}

##
write_boundaries <- function(binaries = binaries){
  
  bounds <- c()
  for(i in 1:length(binaries[[1]])){
    #bounds <- c(bounds, paste("0 <= ", binaries[[1]][i], " <= 1", sep = ""))
    bounds <- c(bounds, paste("\t", binaries[[1]][i], " >= 0\t \t", sep = ""))
    bounds <- c(bounds, paste("\t", binaries[[1]][i], " <= 1\t \t", sep = ""))
  }
  
  return(bounds)
  
}

##
write_equality_constraints <- function(dataMatrix = dataMatrix, binaries = binaries){
  
  dM <- dataMatrix[[1]]
  equalityConstraints <- c()
  if(nrow(dM) > 1){
    for(i in 1:ncol(dM)){
      for(j in 1:(nrow(dM)-1)){
        equalityConstraints <- c(equalityConstraints, paste("xb", i, " - xb", j + ncol(dM) + i - 1, " = 0", sep = ""))
      }
    }
  }
  
  return(equalityConstraints)
  
}

##
write_constraints_1 <- function(dataMatrix = dataMatrix, binaries = binaries, pknList = pknList){
  
  sif <- createSIF(pknList)
  
  constraints1 <- c()
  
  for(i in 1:nrow(sif)){
    
    ss <- which(dataMatrix$species == sif[i, 1])
    tt <- which(dataMatrix$species == sif[i, 3])
    
    yC <- binaries[[1]][dim(dataMatrix[[1]])[1]*dim(dataMatrix[[1]])[2]+i]
    
    constraints1 <- c(constraints1, paste("xb", ss, " + ", "xb", tt, " - 2 ", yC, " >= 0", sep = ""))
    if(sif[i, 1] %in% dataMatrix$species[dataMatrix$tgID]){
      
      constraints1 <- c(constraints1, paste("xb", ss, " + ", "xb", tt, " - ", yC, " >= 0", sep = ""))
      constraints1 <- c(constraints1, paste("xb", ss, " + ", "xb", tt, " - ", yC, " <= 1", sep = ""))
      
    }
    
  }
  
  return(constraints1)
  
}

##
write_constraints_2 <- function(dataMatrix = dataMatrix, binaries = binaries, pknList = pknList){
  
  tNames <- dataMatrix$species[dataMatrix$tgID]
  
  sif <- createSIF(pknList)
  
  tAdjacent <- list()
  for(i in 1:length(tNames)){
    
    temp <- c()
    tAdjacent[[length(tAdjacent)+1]] <- c(temp, which(sif[, 1]==tNames[i]))
    
  }
  
  constraints2 <- c()
  for(i in 1:length(tAdjacent)){
    
    temp <- ""
    
    for(j in 1:length(tAdjacent[[i]])){
      
      if(j==1){
        
        temp <- paste(temp, binaries[[1]][tAdjacent[[i]][j] + dim(dataMatrix[[1]])[1]*dim(dataMatrix[[1]])[2]], sep = "")
        
      }
      else{
        
        temp <- paste(temp, " + ", binaries[[1]][tAdjacent[[i]][j] + dim(dataMatrix[[1]])[1]*dim(dataMatrix[[1]])[2]], sep = "")
        
      }
      
    }
    
    constraints2 <- c(constraints2, paste(temp, " >= 1", sep = ""))
    
  }
  
  for(i in 1:length(tNames)){
    
    constraints2 <- c(constraints2, paste(binaries[[1]][which(dataMatrix$species==tNames[i])], " = 1", sep = ""))
    
  }
  
  return(constraints2)
  
}

##
write_constraints_3 <- function(dataMatrix = dataMatrix, binaries = binaries, pknList = pknList){
  
  nNames <- dataMatrix$species[dataMatrix$dnID]
  
  sif <- createSIF(pknList)
  
  nAdjacent <- list()
  for(i in 1:length(nNames)){
    
    temp <- c()
    nAdjacent[[length(nAdjacent)+1]] <- c(temp, which(sif[, 1]==nNames[i]))
    
  }
  
  constraints3 <- c()
  for(i in 1:length(nAdjacent)){
    temp <- ""
    for(j in 1:length(nAdjacent[[i]])){
      
      if(j==1){
        temp <- paste(temp, binaries[[1]][nAdjacent[[i]][j] + dim(dataMatrix[[1]])[1]*dim(dataMatrix[[1]])[2]], sep = "")
      }
      else{
        temp <- paste(temp, " + ", binaries[[1]][nAdjacent[[i]][j] + dim(dataMatrix[[1]])[1]*dim(dataMatrix[[1]])[2]], sep = "")
      }
      
    }
    
    bb <- which(dataMatrix$species == nNames[i])
    constraints3 <- c(constraints3, paste(temp, " - ", binaries[[1]][bb], " >= 0"))
  }
  
  return(constraints3)
  
}

##
write_constraints_4 <- function(dataMatrix = dataMatrix, binaries = binaries, pknList = pknList){
  
  nNames <- dataMatrix$species[dataMatrix$dnID]
  
  sif <- createSIF(pknList)
  
  nIncident <- list()
  for(i in 1:length(nNames)){
    
    temp <- c()
    nIncident[[length(nIncident)+1]] <- c(temp, which(sif[, 3]==nNames[i]))
    
  }
  
  constraints4 <- c()
  for(i in 1:length(nIncident)){
    if(length(nIncident[[i]]) > 0){
      temp <- ""
      for(j in 1:length(nIncident[[i]])){
        
        if(j==1){
          temp <- paste(temp, binaries[[1]][nIncident[[i]][j] + dim(dataMatrix[[1]])[1]*dim(dataMatrix[[1]])[2]], sep = "")
        }
        else{
          temp <- paste(temp, " + ", binaries[[1]][nIncident[[i]][j] + dim(dataMatrix[[1]])[1]*dim(dataMatrix[[1]])[2]], sep = "")
        }
        
      }
      
      bb <- which(dataMatrix$species == nNames[i])
      constraints4 <- c(constraints4, paste(temp, " - ", binaries[[1]][bb], " >= 0"))
    }
    else{
      
      bb <- which(dataMatrix$species == nNames[i])
      constraints4 <- c(constraints4, paste(binaries[[1]][bb], " = 0"))
      
    }
  }
  
  return(constraints4)
  
}

##
write_constraints_5 <- function(dataMatrix = dataMatrix, binaries = binaries, pknList = pknList){
  
  nNames <- dataMatrix$species[dataMatrix$dsID]
  
  sif <- createSIF(pknList)
  
  nIncident <- list()
  for(i in 1:length(nNames)){
    
    temp <- c()
    nIncident[[length(nIncident)+1]] <- c(temp, which(sif[, 3]==nNames[i]))
    
  }
  
  constraints5 <- c()
  for(i in 1:length(nIncident)){
    if(length(nIncident[[i]]) > 0){
      temp <- ""
      for(j in 1:length(nIncident[[i]])){
        
        if(j==1){
          temp <- paste(temp, binaries[[1]][nIncident[[i]][j] + dim(dataMatrix[[1]])[1]*dim(dataMatrix[[1]])[2]], sep = "")
        }
        else{
          temp <- paste(temp, " + ", binaries[[1]][nIncident[[i]][j] + dim(dataMatrix[[1]])[1]*dim(dataMatrix[[1]])[2]], sep = "")
        }
        
      }
      
      bb <- which(dataMatrix$species == nNames[i])
      constraints5 <- c(constraints5, paste(temp, " - ", binaries[[1]][bb], " >= 0"))
    }
  }
  
  return(constraints5)
  
}


##
all_constraints <- function(equalityConstraints = equalityConstraints, constraints1 = constraints1, constraints2 = constraints2,
                            constraints3 = constraints3, constraints4 = constraints4, constraints5 = constraints5){
  
  kk <- 1
  allConstraints <- c()
  
  if(length(equalityConstraints) > 0){
    
    for(i in 1:length(equalityConstraints)){
      
      allConstraints <- c(allConstraints, paste("c", kk, ":\t", equalityConstraints[i], "\t \t", sep = ""))
      kk <- kk + 1
      
    }
    
  }
  
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
sif2graph<-function(sif){
  
  #if the input is a character it shoud be the name ot the sif file
  #otherwise a matrix in the sif format
  if (is.vector(sif) && (typeof(sif) == "character")){
    sif = read.table(sif) 
  }
  
  # build the unique vertices from the column 1 and 3 of the SIF file
  vertices = unique(c(as.character(sif[,1]), as.character(sif[,3])))
  # some aliases
  v1 = sif[,1]
  v2 = sif[,3]
  edges = as.numeric(sif[,2])
  
  l = length(vertices) - 1
  g <- new("graphNEL", nodes=vertices, edgemode="directed")
  #weights = rep(1, l)
  weights = edges
  for (i in 1:length(v1)){
    g <- addEdge(as.character(v1[i]), as.character(v2[i]), g, weights[i])
  }
  return(g)
  
}

##
readOutSIF<- function(cplexSolutionFileName, binaries = binaries){
  
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
    
    if(round(as.numeric(cplexSolutionEdges[[i]]$variable[3])) == 1){
      
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
reducedModel <- function(sif = sif, dataMatrix = dataMatrix){
  
  gg <- igraph.from.graphNEL(sif2graph(sif = sif))
  
  tNames <- colnames(dataMatrix[[1]])[dataMatrix$tgID]
  sNames <- colnames(dataMatrix[[1]])[dataMatrix$dsID]
  
  allTargets <- c()
  for(i in 1:length(tNames)){
    
    allTargets <- c(allTargets, strsplit(tNames[i], split = ":")[[1]][2])
    
  }
  
  allMeasurements <- c()
  for(i in 1:length(sNames)){
    
    allMeasurements <- c(allMeasurements, strsplit(sNames[i], split = ":")[[1]][2])
    
  }
  
  allSpecies <- rownames(get.adjacency(gg))
  allMeasurements <- allMeasurements[allMeasurements %in% allSpecies]
  
  allPaths <- list()
  for(i in 1:length(allTargets)){
    
    for(j in 1:length(allMeasurements)){
      
      allPaths[[length(allPaths)+1]] <- get.all.shortest.paths(gg, from = which(allSpecies==allTargets[i]), to = which(allSpecies==allMeasurements[j]))
      
    }
  }
  
  sifNew <- matrix(, nrow = 1, ncol = 3)
  sifNew[1, 1] <- allSpecies[allPaths[[1]]$res[[1]][1]]
  sifNew[1, 2] <- 1
  sifNew[1, 3] <- allSpecies[allPaths[[1]]$res[[1]][2]]
  for(i in 1:length(allPaths)){
    
    for(j in 1:length(allPaths[[i]]$res)){
      
      for(k in 1:(length(allPaths[[i]]$res[[j]])-1)){
        
        temp <- matrix(, nrow = 1, ncol = 3)
        temp[1, 1] <- allSpecies[allPaths[[i]]$res[[j]][k]]
        temp[1, 2] <- 1
        temp[1, 3] <- allSpecies[allPaths[[i]]$res[[j]][k+1]]
        
        sifNew <- rbind(sifNew, temp)
        
      }
      
    }
    
  }
  
  return(unique(sifNew))
  
}

##
writeFile <- function(objectiveFunction, constraints, bounds, binaries){
  
  data = "testFile.lp"
  write("enter Problem", data)
  write("", data, append = TRUE)
  write("Minimize", data, append = TRUE)
  
  write(objectiveFunction, data, append = TRUE)
  
  write("Subject To", data, append = TRUE)
  write(constraints, data, append = TRUE)
  
  write("End", data, append = TRUE)
  
}

##
removeRedundancies <- function(sif = sif, dataMatrix = dataMatrix){
  
  tNames <- dataMatrix$species[dataMatrix$tgID]
  sNames <- dataMatrix$species[dataMatrix$dsID]
  
  nodesSIF <- unique(c(sif[, 1], sif[, 3]))
  sNames <- sNames[which(sNames %in% nodesSIF)]
  
  adj <- matrix(0, nrow = length(nodesSIF), ncol = length(nodesSIF))
  rownames(adj) <- nodesSIF
  colnames(adj) <- nodesSIF
  
  for(i in 1:nrow(sif)){
    
    adj[which(rownames(adj) == sif[i, 1]), which(colnames(adj) == sif[i, 3])] <- 1
    
  }
  
  gg <- graph_from_adjacency_matrix(adjmatrix = adj, mode = "directed")
  
  sP <- list()
  for(i in 1:length(tNames)){
    
    for(j in 1:length(sNames)){
      
      allP <- all_simple_paths(gg, from = which(rownames(adj) == tNames[i]), to = which(colnames(adj) == sNames[j]))
      ss <- c()
      for(k in 1:length(allP)){
        
        if(length(allP[[k]]) > 2){
          
          sum = 0
          for(l in 2:(length(allP[[k]])-1)){
            
            currNode <- rownames(adj)[allP[[k]][l]]
            approxNode <- sif[which(sif[, 1]==currNode), 3]
            approxSinkNode <- which(!(approxNode %in% sif[, 1]))
            sum = sum + length(approxSinkNode)
            
          }
          
        }
        
        ss <- c(ss, sum)
        
      }
      
      for(l in 1:length(which(ss==max(ss)))){
        
        sP[[length(sP)+1]] <- unlist(allP[[which(ss==max(ss))[l]]])
        
      }
      
    }
    
  }
  
  allSpecies <- rownames(adj)
  allPaths <- sP
  
  
  sifNew <- matrix(, nrow = 1, ncol = 3)
  sifNew[1, 1] <- allSpecies[allPaths[[1]][1]]
  sifNew[1, 2] <- 1
  sifNew[1, 3] <- allSpecies[allPaths[[1]][2]]
  for(i in 1:length(allPaths)){
    
    for(j in 1:(length(allPaths[[i]])-1)){
      
      temp <- matrix(, nrow = 1, ncol = 3)
      temp[1, 1] <- allSpecies[allPaths[[i]][j]]
      temp[1, 2] <- 1
      temp[1, 3] <- allSpecies[allPaths[[i]][j+1]]
      
      sifNew <- rbind(sifNew, temp)
      
    }
    
  }
  
  return(unique(sifNew))
  
}
