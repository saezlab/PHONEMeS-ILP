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
  
  SIF <- unique(SIF)
  
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

create_binary_variables_for_z_vector <- function(pknList = pknList, dataMatrix = dataMatrix){
  
  sif <- createSIF(pknList = pknList)
  
  z_vector <- paste0(sif[, 1], "=", sif[, 3])
  
  n_interactions <- nrow(sif)
  n_experiments <- nrow(dataMatrix$dataMatrix)
  
  # auxiliary index vectors to use vectorized code instead of nested loops
  index_j <- rep(c(1:n_interactions), n_experiments)
  index_i <- rep(c(1:n_experiments), n_interactions)
  index_i <- c(matrix(index_i, nrow = n_interactions, byrow = TRUE))  # reshape into correct order
  # variable names
  variables <- paste0("z_", index_j, "^", index_i)
  # explanation for each variable
  identifiers <- paste0("interaction ", z_vector, " in experiment ", index_i)
  
  binary_variables_list = list(1:length(variables),variables, identifiers)
  
  return(binary_variables_list)
}

##
create_binary_variables_for_x_vector <- function(dataMatrix = dataMatrix){
  # nrows: number of experiments
  # ncols: number of species
  dim_data <- dim(dataMatrix$dataMatrix)
  n_experiments <- dim_data[1]
  n_species <- dim_data[2]
  
  # names of species
  allSpecies <- gsub("^.*:", "", colnames(dataMatrix$dataMatrix))
  
  # auxiliary index vectors to use vectorized code instead of nested loops
  index_j <- rep(c(1:n_species), n_experiments)
  index_i <- rep(c(1:n_experiments), n_species)
  index_i <- c(matrix(index_i, nrow = n_species, byrow = TRUE))  # reshape into correct order
  # variable names
  variables <- paste0("x_", index_j, "^", index_i)
  # explanation for each variable
  identifiers <- paste("species", allSpecies[index_j], "in experiment", index_i)
  
  binary_variables_list = list(1:length(variables), variables, identifiers)
  
  return(binary_variables_list)
}

##
create_orin_orout_variables <- function(dataMatrix = dataMatrix, targets = targets){
  
  dnID <- c()
  for(ii in 1:nrow(dataMatrix[[1]])){dnID <- c(dnID, dataMatrix$species[dataMatrix$dnID[[ii]]])}
  dnID <- unique(dnID)
  dnID <- setdiff(dnID, unique(unlist(targets)))
  
  variables <- paste0("orin_", 1:length(dnID))
  variables <- c(variables, paste0("orout_", 1:length(dnID)))
  
  identifiers <- paste0("orin ", dnID)
  identifiers <- c(identifiers, paste0("orout ", dnID))
  
  binary_variables_list = list(1:length(variables),variables, identifiers)
  
  return(binary_variables_list)
  
}

##
create_binaries <- function(binaries_x = binaries_x, binaries_z = binaries_z, binaries_in_out = binaries_in_out, binaries_y = binaries_y){
  
  numbers <- c()
  bins <- append(binaries_x[[1]], binaries_z[[1]]+length(binaries_x[[1]]))
  bins <- append(bins, binaries_in_out[[1]]+length(bins))
  bins <- append(bins, binaries_y[[1]]+length(bins))
  numbers <- paste0("xb", bins)
  variables <- append(binaries_x[[2]], binaries_z[[2]])
  variables <- append(variables, binaries_in_out[[2]])
  variables <- append(variables, binaries_y[[2]])
  identifiers <- append(binaries_x[[3]], binaries_z[[3]])
  identifiers <- append(identifiers, binaries_in_out[[3]])
  identifiers <- append(identifiers, binaries_y[[3]])
  
  binaries = list(numbers, variables, identifiers)
  
  return(binaries)
  
}

##
write_objective_function <- function(dataMatrix = dataMatrix, binaries = binaries, sizePen = TRUE, penMode = "edge"){
  
  # # Original set of nested loops
  # dM <- dataMatrix[[1]]
  # kk <- 1
  # objectiveFunction <- "obj:\t"
  # for(i in 1:nrow(dM)){
  #   for(j in 1:ncol(dM)){
  #     if(as.numeric(dM[i, j]) <= 0){
  #       if(as.numeric(dM[i, j]) == 0){
  #         objectiveFunction <- paste(objectiveFunction, "", sep = "")
  #         kk <- kk+1
  #       }
  #       else{
  #         objectiveFunction <- paste(objectiveFunction, " ", 10*as.numeric(dM[i, j]), " ", "xb", kk, sep = "")
  #         kk <- kk+1
  #       }
  #     }
  #     else{
  #       objectiveFunction <- paste(objectiveFunction, " + ", 10*as.numeric(dM[i, j]), " ", "xb", kk, sep = "")
  #       kk <- kk+1
  #     }
  #   }
  # }
  
  # # Alternative 1: setting the format of the numbers with compulsory sign
  # # Doesn't exactly give the same result as the original nested loops because of number formatting
  # dM <- t(dataMatrix[[1]])
  # objectiveFunction <- "obj:\t"
  # index_linear <- which(dM!=0)
  # obj_1 <- paste(sprintf("%+.13f xb%.0f",
  #                        10*as.numeric(dM[index_linear]),
  #                        index_linear),
  #                collapse=" ")
  # objectiveFunction <- paste(objectiveFunction,
  #                            obj_1,
  #                            sep=" ")
  
  # Alternative 2: Fix the consecutive "+" and "-" signs when they appear.
  # Note that we take the transpose of the matrix so that the linear index matches the order in the original nested loops
  dM <- t(dataMatrix$dataMatrix)
  objectiveFunction <- "obj:\t"
  index_linear <- which(dM!=0)
  obj_1 <- paste(" + ",
                 10*as.numeric(dM[index_linear]),
                 " ",
                 "xb", index_linear,
                 sep="",
                 collapse="")
  obj_1 <- gsub("\\+ -", "-", obj_1, perl=FALSE)
  objectiveFunction <- paste(objectiveFunction,
                             obj_1,
                             sep="")
  
  if(sizePen){
    
    if(penMode=="edge"){
      
      # for(i in 1:length(binaries[[3]])){
      #   
      #   if(strsplit(binaries[[3]][i], split = " ")[[1]][1] == "reaction"){
      #     
      #     objectiveFunction <- paste(objectiveFunction, " + 0.001 ", binaries[[1]][i], sep = "")
      #     
      #   }
      #   
      #   if(strsplit(binaries[[3]][i], split = " ")[[1]][1] == "species"){
      #     
      #     objectiveFunction <- paste(objectiveFunction, " + 0.001 ", binaries[[1]][i], sep = "")
      #     
      #   }
      #   
      # }
      
      # The nested loops above can be replaced by this vectorized code
      penalty_reaction <- 0.001
      penalty_species <- 0.001
      is_reaction_var <- grepl("reaction", binaries[[3]])
      is_species_var <- grepl("species", binaries[[3]])
      penalty_vector <- penalty_reaction*is_reaction_var + penalty_species*is_species_var
      is_nonzero <- penalty_vector!=0
      obj_2 <- paste(" + ", penalty_vector[is_nonzero], " ", binaries[[1]][is_nonzero], sep="", collapse="")
      obj_2 <- gsub("\\+ -", "-", obj_2, perl=FALSE)  # in case some of the penalties was negative (which normally wouldn't)
    } else {
      
      # for(i in 1:length(binaries[[3]])){
      #   
      #   if(strsplit(binaries[[3]][i], split = " ")[[1]][1] == "species"){
      #     
      #     objectiveFunction <- paste(objectiveFunction, " + 0.001 ", binaries[[1]][i], sep = "")
      #     
      #   }
      #   
      # }
      
      penalty_species <- 0.001
      is_species_var <- grepl("species", binaries[[3]])
      penalty_vector <- penalty_species*is_species_var
      is_nonzero <- penalty_vector!=0
      obj_2 <- paste(" + ", penalty_vector[is_nonzero], " ", binaries[[1]][is_nonzero], sep="", collapse="")
      obj_2 <- gsub("\\+ -", "-", obj_2, perl=FALSE)  # in case some of the penalties was negative (which normally wouldn't)
    }
    
  } else {
    
    # for(i in 1:length(binaries[[3]])){
    #   
    #   objectiveFunction <- paste(objectiveFunction, " - 0.001 ", binaries[[1]][i], sep = "")
    #   
    # }
    
    penalty_general <- -0.001
    obj_2 <- paste(" + ", penalty_general, " ", binaries[[1]], sep="", collapse="")
    obj_2 <- gsub("\\+ -", "- ", obj_2, perl=FALSE)  # in case some of the penalties was negative
  }
  objectiveFunction <- paste(objectiveFunction,
                             obj_2,
                             sep="")
  
  return(objectiveFunction)
  
}

##
write_boundaries <- function(binaries = binaries, pknList = pknList, M = M, dataMatrix = dataMatrix){
  
  # bounds <- c()
  # for(i in 1:length(binaries[[1]])){
  #   #bounds <- c(bounds, paste("0 <= ", binaries[[1]][i], " <= 1", sep = ""))
  #   bounds <- c(bounds, paste("\t", binaries[[1]][i], " >= 0\t \t", sep = ""))
  #   bounds <- c(bounds, paste("\t", binaries[[1]][i], " <= 1\t \t", sep = ""))
  # }
  # lower bound (lb)
  bounds_lb <- c(paste0("\t", binaries[[1]], " >= 0\t \t"))
  # lower bound (lb)
  bounds_ub <- c(paste0("\t", binaries[[1]], " <= 1\t \t"))
  # merge and reshape to match the output given by the loop
  bounds <- c(t(cbind(bounds_lb, bounds_ub)))
  
  
  sif <- createSIF(pknList = pknList)
  
  species <- unique(c(sif[, 1], sif[, 3]))
  
  # distVar <- c()
  # for(ii in 1:nrow(dataMatrix[[1]])){
  #   distVar <- c(distVar, paste0("dist{", species, "}_", ii))
  # }
  n_experiments <- nrow(dataMatrix$dataMatrix)
  n_distances <- length(species)
  subindex_experiment <- rep(c(1:n_experiments), n_distances)
  subindex_experiment <- c(matrix(subindex_experiment, nrow = n_distances, byrow = TRUE))  # reshape into correct order
  distVar <- paste0("dist{", species, "}_", subindex_experiment)
  
  bounds <- c(bounds, paste0("\t", "0 <= ", distVar, " <= ", M))
  
  return(bounds)
  
}


##
write_equality_constraints <- function(dataMatrix = dataMatrix, binaries = binaries, pknList = pknList){
  
  # dM <- dataMatrix[[1]]
  # equalityConstraints <- c()
  # sif <- createSIF(pknList = pknList)
  # 
  # for(ii in 1:nrow(dM)){
  #   
  #   for(jj in 1:nrow(sif)){
  #     
  #     aa <- binaries[[1]][which(binaries[[3]]==paste0("reaction ", sif[jj, 1], "=", sif[jj, 3]))]
  #     bb <- binaries[[1]][which(binaries[[3]]==paste0("interaction ", sif[jj, 1], "=", sif[jj, 3], " in experiment ", ii))]
  #     
  #     equalityConstraints <- c(equalityConstraints, paste0(aa, " - ", bb, " >= 0"))
  #     
  #   }
  #   
  # }
  
  dM <- dataMatrix$dataMatrix
  sif <- createSIF(pknList = pknList)
  
  n_experiments <- nrow(dM)
  n_interaction <- nrow(sif)
  
  # auxiliary index vectors to use vectorized code instead of nested loops
  index_j <- rep(c(1:n_interaction), n_experiments)
  index_i <- rep(c(1:n_experiments), n_interaction)
  index_i <- c(matrix(index_i, nrow = n_interaction, byrow = TRUE))  # reshape into correct order
  
  description_a <- paste0("reaction ", sif[index_j, 1], "=", sif[index_j, 3])
  description_b <- paste0("interaction ", sif[index_j, 1], "=", sif[index_j, 3], " in experiment ", index_i)
  
  variable_a <- binaries[[1]][match(description_a, binaries[[3]])]
  variable_b <- binaries[[1]][match(description_b, binaries[[3]])]
  
  equalityConstraints1 <- paste0(variable_a, " - ", variable_b, " >= 0")
  
  
  # if(nrow(dM) > 1){
  # 
  #   for(ii in 1:(nrow(dM)-1)){
  # 
  #     for(jj in 1:nrow(sif)){
  # 
  #       aa <- binaries[[1]][which(binaries[[3]]==paste0("interaction ", sif[jj, 1], "=", sif[jj, 3], " in experiment ", ii))]
  #       bb <- binaries[[1]][which(binaries[[3]]==paste0("interaction ", sif[jj, 1], "=", sif[jj, 3], " in experiment ", ii+1))]
  # 
  #       equalityConstraints <- c(equalityConstraints, paste0(aa, " - ", bb, " = 0"))
  # 
  #     }
  # 
  #   }
  # 
  # }
  
  return(equalityConstraints)
  
}

##
write_constraints_1 <- function(dataMatrix = dataMatrix, binaries = binaries, pknList = pknList){
  
  sif <- createSIF(pknList)
  
  constraints1 <- c()
  
  for(ii in 1:nrow(sif)){
    
    ss <- sif[ii, 1]
    tt <- sif[ii, 3]
    
    for(jj in 1:nrow(dataMatrix$dataMatrix)){
      
      aa <- binaries[[1]][which(binaries[[3]]==paste0("species ", ss, " in experiment ", jj))]
      bb <- binaries[[1]][which(binaries[[3]]==paste0("species ", tt, " in experiment ", jj))]
      cc <- binaries[[1]][which(binaries[[3]]==paste0("interaction ", ss, "=", tt, " in experiment ", jj))]
      
      c1 <- paste0(aa," + ", bb, " - 2", cc, " >= 0")
      
      constraints1 <- c(constraints1, c1)
      
    }
    
  }
  
  constraints1 <- unique(constraints1)
  
  return(constraints1)
  
}

##
write_constraints_2 <- function(dataMatrix = dataMatrix, binaries = binaries, pknList = pknList){
  
  constraints2 <- c()
  
  for(ii in 1:nrow(dataMatrix$dataMatrix)){
    
    tNames <- dataMatrix$species[dataMatrix$tgID[[ii]]]
    
    sif <- createSIF(pknList)
    
    tAdjacent <- list()
    for(i in 1:length(tNames)){

      temp <- c()
      tAdjacent[[length(tAdjacent)+1]] <- c(temp, which(sif[, 1]==tNames[i]))

    }

    for(i in 1:length(tAdjacent)){

      temp <- ""

      for(j in 1:length(tAdjacent[[i]])){

        if(j==1){

          temp <- binaries[[1]][which(binaries[[3]]==paste0("interaction ", sif[tAdjacent[[i]][j], 1], "=", sif[tAdjacent[[i]][j], 3], " in experiment ", ii))]

        }
        else{

          temp <- paste0(temp, " + ", binaries[[1]][which(binaries[[3]]==paste0("interaction ", sif[tAdjacent[[i]][j], 1], "=", sif[tAdjacent[[i]][j], 3], " in experiment ", ii))])

        }

      }

      constraints2 <- c(constraints2, paste(temp, " >= 1", sep = ""))

    }
    
    for(i in 1:length(tNames)){
      
      constraints2 <- c(constraints2, paste0(binaries[[1]][which(binaries[[3]]==paste0("species ", tNames[i], " in experiment ", ii))], " = 1"))
      
    }
    
  }
  
  return(constraints2)
  
}

##
write_constraints_3 <- function(dataMatrix = dataMatrix, binaries = binaries, pknList = pknList){
  
  constraints4 <- c()
  
  sif <- createSIF(pknList)
  
  for(ii in 1:nrow(dataMatrix$dataMatrix)){
    
    nNames <- dataMatrix$species[dataMatrix$dnID[[ii]]]
    
    nIncident <- list()
    for(i in 1:length(nNames)){
      
      temp <- c()
      nIncident[[length(nIncident)+1]] <- c(temp, which(sif[, 1]==nNames[i]))
      
    }
    
    for(i in 1:length(nIncident)){
      if(length(nIncident[[i]]) > 0){
        temp <- ""
        for(j in 1:length(nIncident[[i]])){
          
          ss <- sif[nIncident[[i]][j], 1]
          tt <- sif[nIncident[[i]][j], 3]
          
          if(j==1){
            temp <- binaries[[1]][which(binaries[[3]]==paste0("interaction ", ss, "=", tt, " in experiment ", ii))]
          }
          else{
            temp <- paste0(temp, " + ", binaries[[1]][which(binaries[[3]]==paste0("interaction ", ss, "=", tt, " in experiment ", ii))])
          }
          
        }
        
        bb <- which(binaries[[3]]==paste0("species ", nNames[i], " in experiment ", ii))
        constraints4 <- c(constraints4, paste(temp, " - ", binaries[[1]][bb], " >= 0"))
      }
      else{
        
        bb <- which(binaries[[3]]==paste0("species ", nNames[i], " in experiment ", ii))
        constraints4 <- c(constraints4, paste0(binaries[[1]][bb], " = 0"))
        
      }
      
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
write_constraints_4 <- function(dataMatrix = dataMatrix, binaries = binaries, pknList = pknList){
  
  constraints4 <- c()
  
  sif <- createSIF(pknList)
  
  for(ii in 1:nrow(dataMatrix$dataMatrix)){
    
    nNames <- dataMatrix$species[dataMatrix$dnID[[ii]]]
    
    nIncident <- list()
    for(i in 1:length(nNames)){
      
      temp <- c()
      nIncident[[length(nIncident)+1]] <- c(temp, which(sif[, 3]==nNames[i]))
      
    }
    
    for(i in 1:length(nIncident)){
      if(length(nIncident[[i]]) > 0){
        temp <- ""
        for(j in 1:length(nIncident[[i]])){
          
          ss <- sif[nIncident[[i]][j], 1]
          tt <- sif[nIncident[[i]][j], 3]
          
          if(j==1){
            temp <- binaries[[1]][which(binaries[[3]]==paste0("interaction ", ss, "=", tt, " in experiment ", ii))]
          }
          else{
            temp <- paste0(temp, " + ", binaries[[1]][which(binaries[[3]]==paste0("interaction ", ss, "=", tt, " in experiment ", ii))])
          }
          
        }
        
        bb <- which(binaries[[3]]==paste0("species ", nNames[i], " in experiment ", ii))
        constraints4 <- c(constraints4, paste(temp, " - ", binaries[[1]][bb], " >= 0"))
      }
      else{
        
        bb <- which(binaries[[3]]==paste0("species ", nNames[i], " in experiment ", ii))
        constraints4 <- c(constraints4, paste0(binaries[[1]][bb], " = 0"))
        
      }
      
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
  
  constraints4 <- c()
  
  sif <- createSIF(pknList)
  
  for(ii in 1:nrow(dataMatrix$dataMatrix)){
    
    nNames <- dataMatrix$species[dataMatrix$dsID[[ii]]]
    
    nIncident <- list()
    for(i in 1:length(nNames)){
      
      temp <- c()
      nIncident[[length(nIncident)+1]] <- c(temp, which(sif[, 3]==nNames[i]))
      
    }
    
    for(i in 1:length(nIncident)){
      if(length(nIncident[[i]]) > 0){
        temp <- ""
        for(j in 1:length(nIncident[[i]])){
          
          ss <- sif[nIncident[[i]][j], 1]
          tt <- sif[nIncident[[i]][j], 3]
          
          if(j==1){
            temp <- binaries[[1]][which(binaries[[3]]==paste0("interaction ", ss, "=", tt, " in experiment ", ii))]
          }
          else{
            temp <- paste0(temp, " + ", binaries[[1]][which(binaries[[3]]==paste0("interaction ", ss, "=", tt, " in experiment ", ii))])
          }
          
        }
        
        bb <- which(binaries[[3]]==paste0("species ", nNames[i], " in experiment ", ii))
        constraints4 <- c(constraints4, paste(temp, " - ", binaries[[1]][bb], " >= 0"))
      }
      else{
        
        bb <- which(binaries[[3]]==paste0("species ", nNames[i], " in experiment ", ii))
        constraints4 <- c(constraints4, paste0(binaries[[1]][bb], " = 0"))
        
      }
      
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
all_constraints <- function(equalityConstraints = equalityConstraints, constraints1 = constraints1, constraints2 = constraints2,
                            constraints3 = constraints3, constraints4 = constraints4, constraints5 = constraints5, constraints6 = constraints6){
  
  kk <- 1
  allConstraints <- c()
  
  if(length(equalityConstraints)>0){
    
    for(i in 1:length(equalityConstraints)){
      
      allConstraints <- c(allConstraints, paste("c", sprintf("%.0f", kk), ":\t", equalityConstraints[i], "\t \t", sep = ""))
      kk <- kk + 1
      
    }
    
  }
  
  for(i in 1:length(constraints1)){
    
    allConstraints <- c(allConstraints, paste("c", sprintf("%.0f", kk), ":\t", constraints1[i], "\t \t", sep = ""))
    kk <- kk + 1
    
  }
  
  for(i in 1:length(constraints2)){
    
    allConstraints <- c(allConstraints, paste("c", sprintf("%.0f", kk), ":\t", constraints2[i], "\t \t", sep = ""))
    kk <- kk + 1
    
  }
  
  for(i in 1:length(constraints3)){
    
    allConstraints <- c(allConstraints, paste("c", sprintf("%.0f", kk), ":\t", constraints3[i], "\t \t", sep = ""))
    kk <- kk + 1
    
  }
  
  for(i in 1:length(constraints4)){
    
    allConstraints <- c(allConstraints, paste("c", sprintf("%.0f", kk), ":\t", constraints4[i], "\t \t", sep = ""))
    kk <- kk + 1
    
  }
  
  for(i in 1:length(constraints5)){
    
    allConstraints <- c(allConstraints, paste("c", sprintf("%.0f", kk), ":\t", constraints5[i], "\t \t", sep = ""))
    kk <- kk + 1
    
  }
  
  for(i in 1:length(constraints6)){
    
    allConstraints <- c(allConstraints, paste("c", sprintf("%.0f", kk), ":\t", constraints6[i], "\t \t", sep = ""))
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
  
  distVar <- c()
  for(ii in 1:nrow(dataMatrix[[1]])){
    distVar <- c(distVar, paste0("dist{", species, "}_", ii))
  }
  
  speciesVar <- binaries[[1]][grepl(pattern = "species", x = binaries[[3]])]
  speciesExp <- binaries[[3]][grepl(pattern = "species", x = binaries[[3]])]
  
  for(i in 1:length(distVar)){
    
    ss <- speciesVar[grepl(pattern = paste0("species ", species[i], " "), x = speciesExp)]
    
    cc1 <- c(cc1, paste0(ss, " - ", distVar[i], " <= 0"))
    
  }
  
  ##
  reacVar <- binaries[[1]][grepl(pattern = "interaction", x = binaries[[3]])]
  reacExp <- binaries[[3]][grepl(pattern = "interaction", x = binaries[[3]])]
  cnt <- 1
  for(ii in 1:nrow(dataMatrix[[1]])){
    
    for(i in 1:nrow(sif)){
      
      ss <- sif[i, 1]
      tt <- sif[i, 3]
      
      cc2 <- c(cc2, paste0(paste0("dist{", tt, "}_", ii), " - ", paste0("dist{", ss, "}_", ii), " - ", M, " ", reacVar[cnt], " >= ", 1-M))
      
      cnt <- cnt + 1
      
    }
    
  }
  
  ##
  cc3 <- paste0(distVar, " <= ", M)
  
  return(c(cc1, cc2, cc3))
  
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
write_in_out_constraints <- function(binaries = binaries, targets = targets, dataMatrix = dataMatrix, pknList = pknList){
  
  constraints <- c()

  sif <- createSIF(pknList = pknList)

  #
  dnID <- c()
  for(ii in 1:nrow(dataMatrix[[1]])){dnID <- c(dnID, dataMatrix$species[dataMatrix$dnID[[ii]]])}
  dnID <- unique(dnID)
  dnID <- setdiff(dnID, unique(unlist(targets)))

  for(ii in 1:length(dnID)){

    idx1 <- which(sif[, 3]==dnID[ii])
    idx2 <- which(sif[, 1]==dnID[ii])
    orinVar <- binaries[[1]][which(binaries[[3]]==paste0("orin ", dnID[ii]))]
    oroutVar <- binaries[[1]][which(binaries[[3]]==paste0("orout ", dnID[ii]))]

    if((length(idx1)>0) && (length(idx2)>0)){

      varIn <- c()
      for(jj in 1:length(idx1)){

        varIn <- c(varIn, binaries[[1]][which(binaries[[3]]==paste0("reaction ", sif[idx1[jj], 1], "=", sif[idx1[jj], 3]))])

      }
      varOut <- c()
      for(jj in 1:length(idx2)){

        varOut <- c(varOut, binaries[[1]][which(binaries[[3]]==paste0("reaction ", sif[idx2[jj], 1], "=", sif[idx2[jj], 3]))])

      }

      cc1 <- paste0(orinVar, " - ", oroutVar, " = 0")

      cc2 <- ""
      for(kk in 1:length(varIn)){
        if(kk==1){
          cc2 <- paste0(cc2, varIn[kk])
        } else {
          cc2 <- paste0(cc2, " + ", varIn[kk])
        }
      }
      cc2 <- paste0(cc2, " - ", orinVar, " >= 0")

      cc3 <- ""
      for(kk in 1:length(varOut)){
        if(kk==1){
          cc3 <- paste0(cc3, varOut[kk])
        } else {
          cc3 <- paste0(cc3, " + ", varOut[kk])
        }
      }
      cc3 <- paste0(cc3, " - ", oroutVar, " >= 0")

      cc4 <- c()
      for(kk in 1:length(varIn)){
        cc4 <- c(cc4, paste0(varIn[kk], " - ", orinVar, " <= 0"))
      }

      cc5 <- c()
      for(kk in 1:length(varOut)){
        cc5 <- c(cc5, paste0(varOut[kk], " - ", oroutVar, " <= 0"))
      }

      constraints <- c(constraints, c(cc1, cc2, cc3, cc4, cc5))

    }

  }

  return(constraints)
  
}
