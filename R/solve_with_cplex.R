#'\code{solve_with_cplex}
#'
#'Solving problem via CPLEX
#'

solve_with_cplex <- function(solverPath){
  
  if (Sys.info()[1]=="Windows") {
    file.copy(from = solverPath,to = getwd())
    system(paste0("cplex.exe -f cplexCommand_", 
                  condition,"_",repIndex,".txt"))
    file.remove("cplex.exe")
    Elapsed_2 <- proc.time() - ptm
  } else {
    system(paste0(solverPath, " -f cplexCommand_", 
                  condition,"_",repIndex,".txt"))
    Elapsed_2 <- proc.time() - ptm
  }
  
  # load mapping information
  binaries <- readRDS("tmp_binaries.rds")
  
  # Read the results from the CPLEX and do the necessary processing of the model
  library(XML)
  resultsSIF1 <- readOutSIFAll(cplexSolutionFileName = "results1.txt", binaries = binaries)
  colnames(resultsSIF1) <- c("Source", "Interaction", "Target")
  # write.table(resultsSIF1, file = "res1.txt", quote = FALSE, row.names = FALSE, sep = "\t")
  
  # change format to data.frame
  resultsSIF1 <- data.frame(resultsSIF1, stringsAsFactors=FALSE)
  resultsSIF1 <- resultsSIF1 %>% mutate(Interaction=as.numeric(Interaction))
  
  # clean-up temporary files
  file.remove("cplex.log")
  for (i in 0:20)
  {
    file.remove(paste("clone",paste(i,".log", sep = ""), sep = ""))
  }
  file.remove("results1.txt")
  file.remove("testFile.lp")
  file.remove("tmp_binaries.rds")
  file.remove("cplexCommand.txt")
  
  return(resultsSIF1)
}