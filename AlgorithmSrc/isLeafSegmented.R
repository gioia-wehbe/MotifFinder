#GIOIA: My own function to check if node is leaf or not in the case of segmented pstree!!! used in nodeKeepModified
#last updated= 31-01-15

isLeafSegmented <- function(x) 
  {
  test=which(x@leaf==FALSE)
  if(length(test>0))
    return(FALSE)
  else
    return(TRUE)
  }

environment(isLeafSegmented)=asNamespace("PST")