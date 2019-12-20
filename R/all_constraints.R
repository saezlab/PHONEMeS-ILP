#'\code{all_constraints}
#'
#'Grouping all constraints
#'

all_constraints <- function(equalityConstraints = equalityConstraints, 
                            constraints1 = constraints1, 
                            constraints2 = constraints2,
                            constraints3 = constraints3, 
                            constraints4 = constraints4, 
                            constraints5 = constraints5, 
                            constraints6 = constraints6){
  
  kk <- 1
  allConstraints <- c()
  
  if(!is.null(equalityConstraints)){
    if(length(equalityConstraints)>0){
      
      for(i in 1:length(equalityConstraints)){
        
        allConstraints <- 
          c(allConstraints, 
            paste("c", kk, ":\t", equalityConstraints[i], "\t \t", sep = ""))
        kk <- kk + 1
        
      }
      
    }
  }
  
  if(!is.null(constraints1)){
    for(i in 1:length(constraints1)){
      
      allConstraints <- 
        c(allConstraints, 
          paste("c", kk, ":\t", constraints1[i], "\t \t", sep = ""))
      kk <- kk + 1
      
    }
  }
  
  if(!is.null(constraints2)){
    for(i in 1:length(constraints2)){
      
      allConstraints <- 
        c(allConstraints, 
          paste("c", kk, ":\t", constraints2[i], "\t \t", sep = ""))
      kk <- kk + 1
      
    }
  }
  
  if(!is.null(constraints3)){
    for(i in 1:length(constraints3)){
      
      allConstraints <- 
        c(allConstraints, 
          paste("c", kk, ":\t", constraints3[i], "\t \t", sep = ""))
      kk <- kk + 1
      
    }
  }
  
  if(!is.null(constraints4)){
    for(i in 1:length(constraints4)){
      
      allConstraints <- 
        c(allConstraints, 
          paste("c", kk, ":\t", constraints4[i], "\t \t", sep = ""))
      kk <- kk + 1
      
    }
  }
  
  if(!is.null(constraints5)){
    for(i in 1:length(constraints5)){
      
      allConstraints <- 
        c(allConstraints, 
          paste("c", kk, ":\t", constraints5[i], "\t \t", sep = ""))
      kk <- kk + 1
      
    }
  }
  
  if(!is.null(constraints6)){
    for(i in 1:length(constraints6)){
      
      allConstraints <- 
        c(allConstraints, 
          paste("c", kk, ":\t", constraints6[i], "\t \t", sep = ""))
      kk <- kk + 1
      
    }
  }
  
  return(allConstraints)
  
}