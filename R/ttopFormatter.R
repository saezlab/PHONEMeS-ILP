#'\code{ttopFormatter}
#'
#'Preparing ttop on the right output format.
#'
#'Aurelien Dugourd

ttopFormatter <- function(ttop)
{
  ttop$ID <- row.names(ttop)
  ttop <- ttop[,c(7,1,2,3,4,5,6)]
  ttop <- ttop[complete.cases(ttop),]
  return(ttop)
}