#'\code{writeFile}
#'
#'Writing the LP file
#'

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