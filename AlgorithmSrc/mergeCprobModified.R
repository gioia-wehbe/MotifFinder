##################################################################################################
# A function that merges ccounts matrices as needed taking into consideration the difference 
# in row lengths due to the positions of each context that we added
# Source: PST library
# Modified by: Gioia Wehbe
# Last updated: 31-01-15
##################################################################################################

merge.cprob.modified <- function(x,y) 
  {
  if (is.null(x)) 
    {
    res <- y
  } 
  else 
    {
    for (i in names(y)) 
      {
      if (i %in% names(x)) 
        {
        matrix_list=list(x[[i]],y[[i]])
        x[[i]] <- rbind.fill.matrix(matrix_list)
        rownames(x[[i]])=rep(NA,nrow(x[[i]]))
      } 
      else 
        {
        x[[i]] <- y[[i]]
      }
    }
    res <- x
  }
}

environment(merge.cprob.modified)=asNamespace("PST")

