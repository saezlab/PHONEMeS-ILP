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
runPHONEMeS_Downsampling <- function(targets.P, 
                                     conditions, 
                                     inputObj, 
                                     experiments, 
                                     bg, 
                                     nIter = 100, 
                                     nK="all", 
                                     solver="cplex", 
                                     solverPath = "~/Documents/cplex"){
  
  conditions <- conditions[experiments]
  targets <- targets.P[experiments]
  valid_solver_list <- c("cplex", "cbc")
  if (!(solver %in% valid_solver_list)){
    stop(paste0("Select a valid solver option (", paste(valid_solver_list, collapse=", "), ")"))
  }
  
  resList <- list()
  
  for(ii in 1:nIter){
    
    print(paste0("### Iteration ", ii, "/", nIter, "###"))
    
    ss <- sample(x = 1:length(inputObj@res), replace = TRUE)
    ss <- unique(ss)
    temp <- inputObj
    temp@res <- temp@res[ss]
    temp@resFC <- temp@resFC[ss]
    
    data.P <- dataBycond(temp, bg, scaled=TRUE, rowBycond=conditions)
    show(data.P)
    
    speciesP(data.P)
    
    # targets <- targets.P[experiments]
    
    pknList<-build_Nw(data.On=data.P, targets.On=targets, bg=bg, nK=nK)
    
    show(pknList)
    
    TG <- unique(unlist(targets.P))
    
    write_lp_file_1(dataGMM = temp, pknList = pknList, targets = targets, experiments = conditions)
    
    if (solver=="cplex"){
      resultsSIF1 <- solve_with_cplex_tp(solverPath)
    } else if (solver=="cbc"){
      resultsSIF1 <- solve_with_cbc(solverPath)
    } else {
      stop("Select a valid solver option ('cplex', 'cbc')")
    }
    
    resList[[length(resList)+1]] <- resultsSIF1
    
  }
  
  resultsSIF <- combine_networks(resList = resList)
  
  
  # write.table(resultsSIF, file = "resultsSIF.txt", quote = FALSE, row.names = FALSE, sep = "\t")
  return(resultsSIF)
  
}