##################################################################################################
# A function that selects the context with the highest Rel_n out of the redundant contexts.
# it breaks ties randomly if two contexts have the same rel_n. It keeps all possible positions
# of the redundant context as a comma-delimited character string
# Source: PST library
# Modified by: Gioia Wehbe
# Last updated: 31-01-15
##################################################################################################

# library("plyr")
maxRelN<-function(res_matrix_redundant)
{
  all_pos=paste(res_matrix_redundant[,'pos'],collapse=",")
  max=max(res_matrix_redundant[,'rel_n'])
  max_index=which(res_matrix_redundant[,'rel_n']==max)
  if(length(max_index)>0)
  {
    rand=sample(1:length(max_index),size=1)
    max_index=max_index[rand]
  }
  res_matrix_redundant=res_matrix_redundant[max_index,]
  res_matrix_redundant[,'pos']=all_pos
  return(res_matrix_redundant)
}

environment(maxRelN)=asNamespace("PST")


