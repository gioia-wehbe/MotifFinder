#last updated: 31-01-15

getPstrNode<-function(node_name,tree)
{
  #returns the PSTr node with the given node_name from the given PSTf tree
  #NB: it returns NULL if node name is not found or if not in correct format...
  ind=length(unlist(strsplit(node_name,split="-")))
  ind_length_pstrs=tree[[ind+1]]
  pstr_node=ind_length_pstrs[node_name]
  pstr_node=pstr_node[[node_name]]
  return(pstr_node)
}


environment(getPstrNode)=asNamespace("PST")
