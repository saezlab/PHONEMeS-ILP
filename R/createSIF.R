#'\code{createSIF}
#'
#'Creating SIF file from pknList object
#'

createSIF <- function(pknList = pknList){
  
  allInteractions <- pknList@interactions
  SIF <- matrix(, nrow = nrow(allInteractions), ncol = 3)
  
  SIF[, 1] <- allInteractions$K.ID
  
  SIF[, 2] <- rep(1, nrow(allInteractions))
  
  SIF[, 3] <- allInteractions$S.cc
  
  SIF <- unique(SIF)
  
  return(SIF)
  
}