getNodestoKeep<-function(all_nodes_names,tree)
{
  keep_list=c()
  end_index=length(all_nodes_names)
  for(i in 1:end_index)
  {
    current_context=all_nodes_names[i]
    start_index=i+1
    search_space=all_nodes_names[start_index:end_index]
    hits=grep(current_context,search_space)
    if(length(hits)==0)  #if there are no hits => not larger superfix => this is the longest suffix => add to the keep list
    {
      keep_list=c(keep_list,current_context)
    }
  }
  
  keep_list=c(keep_list,current_context)
  keep_list=keep_list[-1]   #remove "e"
  return(keep_list)
}

environment(getNodestoKeep)=asNamespace("PST")
