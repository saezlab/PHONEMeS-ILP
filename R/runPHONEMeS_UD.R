#' Run PHONEMeS ILP
#' 
#' @param Arguments 
#' @param targets.P
#' @param conditions
#' @param dataGMM
#' @param experiments
#' @param bg
#' @param nK 
#' @param solver Solver to use for solving the ILP.
#
#' @return SIF like data.frame with the output network.
runPHONEMeS_UD <- function(targets.P, 
                           conditions, 
                           dataGMM, 
                           experiments, 
                           bg, 
                           nK="all", 
                           solver="cplex", 
                           solverPath="~/Documents/cplex"){
  
  conditions <- conditions[experiments]
  valid_solver_list <- c("cplex", "cbc")
  if (!(solver %in% valid_solver_list)){
    stop(paste0("Select a valid solver option (", paste(valid_solver_list, collapse=", "), ")"))
  }
  
  data.P <- dataBycond(dataGMM, bg, scaled=TRUE, rowBycond=conditions)
  show(data.P)
  
  speciesP(data.P)
  
  targets.P <- targets.P[experiments]
  
  pknList<-build_Nw(data.On=data.P, targets.On=targets.P, bg=bg, nK=nK)
  
  show(pknList)
  
  TG <- unique(unlist(targets.P))
  
  write_lp_file_inv_1(dataGMM = dataGMM, pknList = pknList, targets = targets.P, experiments = conditions)
  
  if (solver=="cplex"){
    resultsSIF1 <- solve_with_cplex(solverPath)
  } else if (solver=="cbc"){
    resultsSIF1 <- solve_with_cbc(solverPath)
  } else {
    stop("Select a valid solver option ('cplex', 'cbc')")
  }
  
  # write.table(resultsSIF, file = "resultsSIF.txt", quote = FALSE, row.names = FALSE, sep = "\t")
  colnames(resultsSIF1) <- c("Source", "Interaction", "Target")
  resultsSIF1[, 2] <- "1"
  
  #####
  # Step - 2
  #Create the PKN list that will be used for optimisation
  pknList <- build_Nw_Inv(data.On=data.P, targets.On=targets.P, bg=bg,nK=nK)
  
  if(!is.null(pknList)){
    show(pknList)
    
    TG <- unlist(targets.P)
    
    write_lp_file_inv_2(dataGMM = dataGMM, pknList = pknList, targets = targets.P, experiments = conditions)
    
    if (solver=="cplex"){
      resultsSIF2 <- solve_with_cplex(solverPath)
      resultsSIF2 <- resultsSIF2[,c(3,2,1)]
    } else if (solver=="cbc"){
      resultsSIF2 <- solve_with_cbc(solverPath)
      resultsSIF2 <- resultsSIF2[,c(3,2,1)]
    } else {
      stop("Select a valid solver option ('cplex', 'cbc')")
    }
    
    # write.table(resultsSIF, file = "resultsSIF.txt", quote = FALSE, row.names = FALSE, sep = "\t")
    colnames(resultsSIF2) <- c("Source", "Interaction", "Target")
    # colnames(resultsSIF2) <- c("Target", "Interaction", "Source")
    resultsSIF2[, 2] <- "1"
    
    resultSIF <- resultsSIF1
    for(ii in 1:nrow(resultsSIF2)){
      
      ss <- resultsSIF2[ii, 1]
      tt <- resultsSIF2[ii, 3]
      idx1 <- which(resultSIF[, 1]==ss)
      idx2 <- which(resultSIF[, 3]==tt)
      idx <- intersect(x = idx1, y = idx2)
      
      if(length(idx)>0){
        resultSIF[idx, 2] <- mean(resultSIF[idx, 2], resultsSIF2[ii, 2])
      } else {
        resultSIF <- rbind(resultSIF, resultsSIF2[ii,])
      }
      
    }
    
    resList <- list()
    resList[[1]] <- resultsSIF1
    # resList[[2]] <- resultsSIF2
    # resList[[3]] <- resultSIF
    resList[[2]] <- resultsSIF2
    # cnt = nrow(resultsSIF2)
    # for(ii in 1:length(cnt)){
    #   resultSIF[nrow(resultsSIF1)+ii, 1] = resultSIF[nrow(resultsSIF1)+ii, 3]
    #   resultSIF[nrow(resultsSIF1)+ii, 3] = resultSIF[nrow(resultsSIF1)+ii, 1]
    # }
    resList[[3]] <- resultSIF
    
    names(resList) <- c("Downside", "Upside", "Combined")
    
    return(resList)
    
  } else {
    
    resList <- list()
    resList[[1]] <- resultsSIF1
    resList[[2]] <- NULL
    resultSIF <- resultsSIF1
    resList[[3]] <- resultSIF
    
    return(resList)
    
  }
  
  # solve_with_cplex <- function(){
  #   system(paste0(getwd(), "/cplex -f cplexCommand.txt"))
  #   
  #   # load mapping information
  #   binaries <- readRDS("tmp_binaries.rds")
  #   
  #   # Read the results from the CPLEX and do the necessary processing of the model
  #   library(XML)
  #   resultsSIF1 <- readOutSIF(cplexSolutionFileName = "results1.txt", binaries = binaries)
  #   colnames(resultsSIF1) <- c("Source", "Interaction", "Target")
  #   # write.table(resultsSIF1, file = "res1.txt", quote = FALSE, row.names = FALSE, sep = "\t")
  #   
  #   # change format to data.frame
  #   resultsSIF1 <- data.frame(resultsSIF1, stringsAsFactors=FALSE)
  #   resultsSIF1 <- resultsSIF1 %>% mutate(Interaction=as.numeric(Interaction))
  #   
  #   # clean-up temporary files
  #   file.remove("cplex.log")
  #   for (i in 0:20)
  #   {
  #     file.remove(paste("clone",paste(i,".log", sep = ""), sep = ""))
  #   }
  #   file.remove("results1.txt")
  #   file.remove("testFile.lp")
  #   file.remove("tmp_binaries.rds")
  #   file.remove("cplexCommand.txt")
  #   
  #   return(resultsSIF1)
  # }
  # 
  # 
  # solve_with_cbc <- function(){
  #   cbc_command <- "cbc testFile.lp solve printi csv solu results_cbc.txt"
  #   system(cbc_command)
  #   
  #   # retrieve solution
  #   readCbcSolution <- function(file, binaries){
  #     library("dplyr")
  #     library("tidyr")
  #     
  #     # read cbc solution file
  #     cbc_table <- read.csv("results_cbc.txt")
  #     
  #     # load mapping information
  #     binaries <- readRDS("tmp_binaries.rds")
  #     
  #     # find mapping of variables from the MILP to the model
  #     mapping_table <- data.frame(binaries)
  #     colnames(mapping_table) <- c("milp", "math", "description")
  #     # keep only reaction variables
  #     mapping_table <- mapping_table %>% filter(grepl("reaction", description))
  #     # merge with solution
  #     cbc_table <- merge(cbc_table, mapping_table, by.x="name", by.y="milp")
  #     # keep only active reactions
  #     cbc_table <- cbc_table %>% mutate(solution=round(as.numeric(solution))) %>% 
  #       filter(solution==1)
  #     # create SIF table
  #     sif <- cbc_table %>% select(description) %>% 
  #       mutate(description = gsub("reaction ", "", description)) %>% 
  #       separate(description, into=c("Source", "Target"), sep="=") %>% 
  #       mutate(Interaction = 1) %>% 
  #       select(Source, Interaction, Target)
  #     
  #     return(sif)
  #   }
  #   
  #   resultsSIF1 <- readCbcSolution("results_cbc.txt", binaries)
  #   
  #   # clean-up temporary files
  #   file.remove("results_cbc.txt")
  #   file.remove("testFile.lp")
  #   file.remove("tmp_binaries.rds")
  #   
  #   return(resultsSIF1)
  # }
  
}
