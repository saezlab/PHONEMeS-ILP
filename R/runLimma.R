#' Running the Limma analysis
#' 
#' Performing differential expression analysis with Limma
#' 
#' @param measurements Phosphoproteomics measurements
#' @param targets Targets of the samples
#' @param comparisons Comparisons for the differential analysis
#' @param pool
#' @param regress_out
#
#' @return Limma output.
#'
#' Aurelien Dugourd
#' 

runLimma <- function(measurements, targets, comparisons = NULL, pool = NULL, regress_out = NULL)
{
  input_check <- checkInputs(measurements, targets)
  if (input_check[[1]]) #input has correct format
  {
    if (!is.null(comparisons))
    {
      if (!is.null(regress_out))
      {
        for (regressor in regress_out)
        {
          measurements <- removeBatchEffect(measurements, targets[,regressor])
        }
      }
      
      cont.matrix <- makeContrastsAlt(targets, comparisons)
      
      if (!is.null(pool))
      {
        cont.matrix <- poolContrasts(cont.matrix, pool)
      }
      
      cont.matrix <- as.data.frame(cont.matrix)
      row.names(cont.matrix) <- unique(targets$condition)
      cont.matrix <- as.matrix(cont.matrix)
      
      fcond <- factor(targets$condition, levels = unique(targets$condition))
      
      design <- model.matrix(~0+fcond)
      design <- as.data.frame(design)
      names(design) <- unique(targets$condition)
      design <- as.matrix(design)
      
      print(cont.matrix)
      
      fit <- lmFit(measurements, design)
      fit2 <- contrasts.fit(fit, cont.matrix)
      fit2 <- eBayes(fit2)
      
      return(list(fit2, cont.matrix, fit))
    }
  }
  else
  {
    print(input_check[[2]])
    return(input_check[[1]])
  }
}