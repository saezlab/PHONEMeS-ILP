#' Run time-point variant PHONEMeS ILP
#' 
#' Main running function for time-point analysis with PHONEMeS.
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
#' @return SIF like data.frame with the combined output network.
#' 

runPHONEMeS_mult <- function(targets.P, conditions, inputObj, experiments, bg, 
                             nIter = 100, nK="all", solver="cplex", 
                             solverPath="~/Documents/cplex"){
  
  valid_solver_list <- c("cplex", "cbc")
  if (!(solver %in% valid_solver_list)){
    stop(paste0("Select a valid solver option (", 
                paste(valid_solver_list, collapse=", "), ")"))
  }
  
  resList <- list()
  
  for(ii in 1:nIter){
    
    print(paste0("### Iteration ", ii, "/", nIter, "###"))
    
    ss <- sample(x = 1:length(inputObj@res), replace = TRUE)
    ss <- unique(ss)
    temp <- inputObj
    temp@res <- temp@res[ss]
    temp@resFC <- temp@resFC[ss]
    
    resListSep <- list()
    
    for(jj in 1:length(experiments)){
      
      if(jj==1){
        
        targets <- targets.P[experiments[[jj]]]
        
        data.P <- dataBycond(temp, bg, scaled=TRUE, 
                             rowBycond=conditions[experiments[[jj]]])
        show(data.P)
        
        speciesP(data.P)
        
        pknList<-build_Nw(data.On=data.P, targets.On=targets, bg=bg, nK=nK)
        pknListTemp <- pknList
        
        show(pknList)
        
        TG <- unique(unlist(targets.P))
        
        write_lp_file_1(dataGMM = temp, pknList = pknList, targets = targets, 
                        experiments = conditions[experiments[[jj]]])
        
        if (solver=="cplex"){
          resultsSIF1 <- solve_with_cplex_tp(solverPath)
        } else if (solver=="cbc"){
          resultsSIF1 <- solve_with_cbc(solverPath)
        } else {
          stop("Select a valid solver option ('cplex', 'cbc')")
        }
        
        resListSep[[length(resListSep)+1]] <- 
          removeRedundantNodes(resultsSIF1 = resultsSIF1)
        
      } else {
        
        targets <- targets.P[experiments[[jj]]]
        
        data.P <- dataBycond(temp, bg, scaled=TRUE, 
                             rowBycond=conditions[experiments[[jj]]])
        show(data.P)
        
        speciesP(data.P)
        
        pknList<-
          build_Nw(data.On=data.P, targets.On=targets.P, bg=bg, nK = "yes")
        pknListTemp@interactions <- 
          unique(rbind(pknList@interactions, pknListTemp@interactions))
        pknListTemp@species <- 
          unique(c(pknList@species, pknListTemp@species))
        
        show(pknListTemp)
        
        write_lp_file_2(prevSIF = resultsSIF1, dataGMM = temp, 
                        pknList = pknListTemp, targets = targets, 
                        experiments = conditions[experiments[[jj]]])
        
        if (solver=="cplex"){
          resultsSIF1 <- solve_with_cplex_tp(solverPath)
        } else if (solver=="cbc"){
          resultsSIF1 <- solve_with_cbc(solverPath)
        } else {
          stop("Select a valid solver option ('cplex', 'cbc')")
        }
        
        
        resListSep[[length(resListSep)+1]] <- 
          removeRedundantNodes(resultsSIF1 = resultsSIF1)
        
      }
      
    }
    
    resList[[length(resList)+1]] <- resListSep
    
  }
  
  sif <- assign_tp_attributes(sifList = resList)
  
  return(sif)
  
}