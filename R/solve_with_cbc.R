#'\code{solve_with_cbc}
#'
#'Solving problem via CPLEX
#'

solve_with_cbc <- function(){
  cbc_command <- paste0(solverPath, 
                        " testFile.lp solve printi csv solu results_cbc.txt")
  system(cbc_command)
  
  # retrieve solution
  readCbcSolution <- function(file, binaries){
    library("dplyr")
    library("tidyr")
    
    # read cbc solution file
    cbc_table <- read.csv("results_cbc.txt")
    
    # load mapping information
    binaries <- readRDS("tmp_binaries.rds")
    
    # find mapping of variables from the MILP to the model
    mapping_table <- data.frame(binaries)
    colnames(mapping_table) <- c("milp", "math", "description")
    # keep only reaction variables
    mapping_table <- mapping_table %>% filter(grepl("reaction", description))
    # merge with solution
    cbc_table <- merge(cbc_table, mapping_table, by.x="name", by.y="milp")
    # keep only active reactions
    cbc_table <- cbc_table %>% mutate(solution=round(as.numeric(solution))) %>% 
      filter(solution==1)
    # create SIF table
    sif <- cbc_table %>% select(description) %>% 
      mutate(description = gsub("reaction ", "", description)) %>% 
      separate(description, into=c("Source", "Target"), sep="=") %>% 
      mutate(Interaction = 1) %>% 
      select(Source, Interaction, Target)
    
    return(sif)
  }
  
  resultsSIF1 <- readCbcSolution("results_cbc.txt", binaries)
  
  # clean-up temporary files
  file.remove("results_cbc.txt")
  file.remove("testFile.lp")
  file.remove("tmp_binaries.rds")
  
  return(resultsSIF1)
}