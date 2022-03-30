## If node not in keep.list it is tagged as pruned
#last updated: 31-01-15

nodeKeepModified <- function(x, keep.list, clist) 
  {
  print("inside nodeKeepModified")
  if (!is.null(clist)) 
    {
#     for(d in 1:length(clist))
#     {
#       print("d")
#       print(d)
#       curr_node=clist[[d]]
#     }
    tmp <- unlist(lapply(clist, node.parent))
  }	
  else
  {
    print("clist is NUUUUUUUUUUULL")
  }
  if (!x@path %in% keep.list && ((!is.null(clist) && isLeafSegmented(x)) | (!is.null(clist) && !x@path %in% tmp))) #GIOIA: did not consider segmented case to check if node is leaf or not!!!
  {
    print("we have pruned current node!!")
    x@pruned[] <- TRUE
  }
else
{
  print("we have NOT pruned current node!!")
}
  return(x)
}

environment(nodeKeepModified)=asNamespace("PST")

