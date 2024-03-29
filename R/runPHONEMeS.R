#' Run PHONEMeS ILP
#' 
#' Main running function for PHONEMeS.
#' 
#' @param Arguments 
#' @param targets.P
#' @param conditions
#' @param dataGMM
#' @param experiments
#' @param bg
#' @param nK 
#' @param solver Solver to use for solving the ILP.
#' @param penFac Penalty factor (default: 0.0001)
#
#' @return SIF like data.frame with the combined output network.
#' 

runPHONEMeS <- function(targets.P, 
                        conditions, 
                        inputObj, 
                        experiments, 
                        bg, 
                        nK="all", 
                        solver="cplex",
                        solverPath = "/usr/bin/cplex",
                        nSolutions=100, 
                        mipgap=0, 
                        relgap=0, 
                        replace=2, 
                        populate=5000, 
                        intensity=4, 
                        timelimit=3600, 
                        penFac = 0.0001){
  
  conditions <- conditions[experiments]
  valid_solver_list <- c("cplex", "cbc")
  if (!(solver %in% valid_solver_list)){
    stop(paste0("Select a valid solver option (", 
                paste(valid_solver_list, collapse=", "), ")"))
  }
  
  data.P <- dataBycond(dataGMM=inputObj, bg, scaled=TRUE, rowBycond=conditions)
  show(data.P)
  
  speciesP(data.P)
  
  targets.P <- targets.P[experiments]
  
  pknList<-build_Nw(data.On=data.P, targets.On=targets.P, bg=bg, nK=nK)
  
  show(pknList)
  
  TG <- unique(unlist(targets.P))
  
  write_lp_file(dataGMM = inputObj, pknList = pknList, targets = targets.P, 
                experiments = conditions, nSolutions=nSolutions, mipgap=mipgap, 
                relgap=relgap, replace=replace, populate=populate, 
                intensity=intensity, timelimit=timelimit, penFac = penFac)
  
  if (solver=="cplex"){
    resultsSIF1 <- solve_with_cplex(solverPath)
  } else if (solver=="cbc"){
    resultsSIF1 <- solve_with_cbc(solverPath)
  } else {
    stop("Select a valid solver option ('cplex', 'cbc')")
  }
  
  colnames(resultsSIF1) <- c("Source", "f50", "Target")
  
  return(resultsSIF1)
  
}