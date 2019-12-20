#'\code{solve_with_cplex_tp}
#'
#'Solving problem via CPLEX - time-point analysis
#'

solve_with_cplex_tp <- function(){
  system(paste0(getwd(), "/cplex -f cplexCommand.txt"))
  
  # load mapping information
  binaries <- readRDS("tmp_binaries.rds")
  
  # Read the results from the CPLEX and do the necessary processing of the model
  library(XML)
  resultsSIF1 <- readOutSIF(cplexSolutionFileName = "results1.txt", 
                            binaries = binaries)
  colnames(resultsSIF1) <- c("Source", "f50", "Target")
  
  # change format to data.frame
  resultsSIF1 <- data.frame(resultsSIF1, stringsAsFactors=FALSE)
  resultsSIF1 <- resultsSIF1 %>% mutate(f50=as.numeric(f50))
  
  # clean-up temporary files
  file.remove("cplex.log")
  for (i in 0:20)
  {
    file.remove(paste("clone",paste(i,".log", sep = ""), sep = ""))
  }
  file.remove("results1.txt")
  file.remove("testFile.lp")
  file.remove("tmp_binaries.rds")
  
  return(resultsSIF1)
}