#
#  This file is part of the CNO software
#
#  Copyright (c) 2018 - RWTH Aachen - JRC COMBINE
#
#  File author(s): E. Gjerga (enio.gjerga@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  CNO website: https://saezlab.github.io/PHONEMeS/
#
##############################################################################
# $Id$

# PHONEMeS - ILP implementation


createSIF <- function(pknList = pknList){
  
  # if((class(pknList)=="PKNlist")==FALSE){
  #   stop("Input object not of class PKNlist")
  # }
  
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
  
  y_vector <- paste0(sif[, 1], "=", sif[, 3])
  
  variables <- paste0("y_", 1:length(y_vector))
  
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
write_objective_function <- function(dataMatrix = dataMatrix, binaries = binaries, sizePen = TRUE, penMode = "edge"){
  
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
          objectiveFunction <- paste(objectiveFunction, " ", 10*as.numeric(dM[i, j]), " ", "xb", kk, sep = "")
          kk <- kk+1
        }
      }
      else{
        objectiveFunction <- paste(objectiveFunction, " + ", 10*as.numeric(dM[i, j]), " ", "xb", kk, sep = "")
        kk <- kk+1
      }
    }
  }
  
  if(sizePen){
    
    if(penMode=="edge"){
      
      for(i in 1:length(binaries[[3]])){
        
        if(strsplit(binaries[[3]][i], split = " ")[[1]][1] == "reaction"){
          
          objectiveFunction <- paste(objectiveFunction, " + 0.001 ", binaries[[1]][i], sep = "")
          
        }
        
      }
      
    } else {
      
      for(i in 1:length(binaries[[3]])){
        
        if(strsplit(binaries[[3]][i], split = " ")[[1]][1] == "species"){
          
          objectiveFunction <- paste(objectiveFunction, " + 0.001 ", binaries[[1]][i], sep = "")
          
        } else {
          
          objectiveFunction <- paste(objectiveFunction, " - 0.0001 ", binaries[[1]][i], sep = "")
          
        }
        
      }
      
    }
    
  } else {
    
    for(i in 1:length(binaries[[3]])){
      
      objectiveFunction <- paste(objectiveFunction, " - 0.001 ", binaries[[1]][i], sep = "")
      
    }
    
  }
  
  return(objectiveFunction)
  
}

##
write_boundaries <- function(binaries = binaries, pknList = pknList, M = M){
  
  bounds <- c()
  for(i in 1:length(binaries[[1]])){
    #bounds <- c(bounds, paste("0 <= ", binaries[[1]][i], " <= 1", sep = ""))
    bounds <- c(bounds, paste("\t", binaries[[1]][i], " >= 0\t \t", sep = ""))
    bounds <- c(bounds, paste("\t", binaries[[1]][i], " <= 1\t \t", sep = ""))
  }
  
  sif <- createSIF(pknList = pknList)
  
  species <- unique(c(sif[, 1], sif[, 3]))
  
  distVar <- paste0("dist{", species, "}")
  
  bounds <- c(bounds, paste0("\t", "0 <= ", distVar, " <= ", M))
  
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
  
  if(length(grep(pattern = "NA", x = constraints3))==0){
    
    return(constraints3)
    
  }
  else{
    
    return(constraints3[-grep(pattern = "NA", x = constraints3)])
    
  }
  
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
  
  if(length(grep(pattern = "NA", x = constraints4))==0){
    
    return(constraints4)
    
  }
  else{
    
    return(constraints4[-grep(pattern = "NA", x = constraints4)])
    
  }
  
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
                            constraints3 = constraints3, constraints4 = constraints4, constraints5 = constraints5, constraints6 = constraints6){
  
  kk <- 1
  allConstraints <- c()
  
  for(i in 1:length(equalityConstraints)){
    
    allConstraints <- c(allConstraints, paste("c", kk, ":\t", equalityConstraints[i], "\t \t", sep = ""))
    kk <- kk + 1
    
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
  
  for(i in 1:length(constraints6)){
    
    allConstraints <- c(allConstraints, paste("c", kk, ":\t", constraints6[i], "\t \t", sep = ""))
    kk <- kk + 1
    
  }
  
  return(allConstraints[3:length(allConstraints)])
  
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
readOutSIFAll<- function(cplexSolutionFileName, binaries = binaries){
  
  reacIDX <- c()
  for(i in 1:length(binaries[[3]])){
    if(strsplit(binaries[[3]][i], split = " ")[[1]][1] == "reaction"){
      reacIDX <- c(reacIDX, i)
    }
  }
  reacVar <- binaries[[1]][reacIDX]
  
  cplexSolutionData <- xmlParse(cplexSolutionFileName)
  cplexSolution <- xmlToList(cplexSolutionData)
  
  cplexSolutionEdgesAll <- list()
  sifAll <- list()
  
  for(ii in 1:(length(cplexSolution)-1)){
    
    sif <- matrix(data = NA, nrow = 1, ncol = 3)
    
    currSolution <- cplexSolution[[ii]][[4]]
    
    for(jj in 1:length(currSolution)){
      
      if((currSolution[[jj]][1]%in%reacVar) && round(as.numeric(currSolution[[jj]][3])==1)){
        
        reaction <- strsplit(binaries[[3]][which(binaries[[1]]==currSolution[[jj]][1])], split = " ")[[1]][2]
        
        sif2bind <- as.matrix(c(strsplit(x = reaction, split = "=")[[1]][1], "1", strsplit(x = reaction, split = "=")[[1]][2]))
        sif2bind <- t(sif2bind)
        sif <- rbind(sif, sif2bind)
        
      }
      
    }
    
    sifAll[[length(sifAll)+1]] <- unique(sif[-1, ])
    
  }
  
  for(ii in 1:length(sifAll)){
    
    if(ii==1){
      
      sif <- sifAll[[1]]
      
    } else {
      
      # sif <- unique(rbind(sif, sifAll[[ii]]))
      for(jj in 1:nrow(sifAll[[ii]])){
        
        idx1 <- which(sif[, 1]==sifAll[[ii]][jj, 1])
        idx2 <- which(sif[, 3]==sifAll[[ii]][jj, 3])
        
        idx <- intersect(x = idx1, y = idx2)
        
        if(length(idx) > 0){
          
          sif[idx, 2] <- as.character(as.numeric(sif[idx, 2])+1)
          
        } else {
          
          sif <- rbind(sif, sifAll[[ii]][jj, ])
          
        }
        
      }
      
    }
    
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

removeRedundantEdges <- function(resultsSIF1 = resultsSIF1){
  
  kinases <- c()
  for(i in 1:nrow(resultsSIF1)){
    
    if(length(strsplit(resultsSIF1[i, 3], split = "_")[[1]])==2){
      
      kinases <- c(kinases, resultsSIF1[i, 3])
      
    }
    
  }
  kinases <- unique(kinases)
  
  toRem <- c()
  for(i in 1:length(kinases)){
    
    idx <- which(resultsSIF1[, 3]==kinases[i])
    
    if(length(idx) > 1){
      
      ctrl <- 0
      for(j in 1:length(idx)){
        
        if(strsplit(resultsSIF1[idx[j], 1], split = "_")[[1]][3] != "R1"){
          
          ctrl <- ctrl +1
          
        }
        
      }
      
      if(ctrl > 0){
        
        for(j in 1:length(idx)){
          
          if(strsplit(resultsSIF1[idx[j], 1], split = "_")[[1]][3] == "R1"){
            
            toRem <- c(toRem, idx[j])
            
          }
          
        }
        
      }
      
    }
    
  }
  
  if(length(toRem) > 0){
    
    idxRem <- c()
    
    speciesRem <- resultsSIF1[toRem, 1]
    
    for(i in 1:nrow(resultsSIF1)){
      
      if(resultsSIF1[i, 1]%in%speciesRem || resultsSIF1[i, 3]%in%speciesRem){
        
        idxRem <- c(idxRem, i)
        
      }
      
    }
    
    resultsSIF1 <- resultsSIF1[-idxRem, ]
    
    return(resultsSIF1)
    
  }
  else{
    
    return(resultsSIF1)
    
  }
  
}

##
removeRedundantNodes <- function(resultsSIF1 = resultsSIF1){
  
  idxToRem <- c()
  
  temp <- resultsSIF1
  
  for(i in 1:nrow(resultsSIF1)){
    
    temp[i, 1] <- gsub(pattern = "_R1", replacement = "", x = resultsSIF1[i, 1])
    temp[i, 3] <- gsub(pattern = "_R1", replacement = "", x = resultsSIF1[i, 3])
    
  }
  
  for(i in 1:nrow(resultsSIF1)){
    
    if(temp[i, 1]==temp[i, 3]){
      
      idxToRem <- c(idxToRem, i)
      
    }
    
  }
  
  return(temp[-idxToRem, ])
  
}

##
write_self_activating_constraints <- function(pknList = pknList, binaries = binaries, dataMatrix = dataMatrix, M = M){
  
  cc1 <- c()
  cc2 <- c()
  cc3 <- c()
  
  sif <- createSIF(pknList = pknList)
  
  species <- unique(c(sif[, 1], sif[, 3]))
  
  distVar <- paste0("dist{", species, "}")
  
  speciesVar <- binaries[[1]][grepl(pattern = "species", x = binaries[[3]])]
  speciesExp <- binaries[[3]][grepl(pattern = "species", x = binaries[[3]])]
  
  for(i in 1:length(distVar)){
    
    ss <- speciesVar[grepl(pattern = paste0("species ", species[i], " "), x = speciesExp)]
    
    cc1 <- c(cc1, paste0(ss, " - ", distVar[i], " <= 0"))
    
  }
  
  ##
  reacVar <- binaries[[1]][grepl(pattern = "reaction", x = binaries[[3]])]
  reacExp <- binaries[[3]][grepl(pattern = "reaction", x = binaries[[3]])]
  for(i in 1:nrow(sif)){
    
    ss <- sif[i, 1]
    tt <- sif[i, 3]
    
    cc2 <- c(cc2, paste0(paste0("dist{", tt, "}"), " - ", paste0("dist{", ss, "}"), " - ", M, " ", reacVar[i], " >= ", 1-M))
    
  }
  
  ##
  cc3 <- paste0(distVar, " <= ", M)
  
  return(c(cc1, cc2, cc3))
  
}