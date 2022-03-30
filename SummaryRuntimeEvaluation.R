#########################################
# Plotting runtime graph
#########################################


print("Plotting Runtime Graph")

param_list=getParamsForPlotTitle()

data_type_dir=""
if(data_type=="Training")
{
  data_type_dir=getTrainingAllOutputDir()
}else if(data_type=="Testing"){
  data_type_dir=getTestingAllOutputDir()
}else {
  data_type_dir=getCompleteDataRuntimeEvaluationDir()
}

#read table of all runs
all_runs_runtime_summary_file=paste(data_type_dir,"/Table_runtime_summary_all_runs.csv",sep="")
all_runs_runtime_summary_table=read.csv(all_runs_runtime_summary_file,sep=",",header=FALSE, colClasses="character",stringsAsFactors=FALSE)
col_names=all_runs_runtime_summary_table[1,]
colnames(all_runs_runtime_summary_table)=col_names
all_runs_runtime_summary_table=all_runs_runtime_summary_table[-1,]

num_code_sec=13
num_measures=4 #Avg,max,min,stdv

# sorting_data_runtime_indexes=seq(from=1,to=nrow(all_runs_runtime_summary_table),by=num_code_sec)
# generating_freq_table_runtime_indexes=seq(from=2,to=nrow(all_runs_runtime_summary_table),by=num_code_sec)
# generating_rel_freq_runtime_indexes=seq(from=3,to=nrow(all_runs_runtime_summary_table),by=num_code_sec)
generating_pvalues_runtime_indexes=seq(from=1,to=nrow(all_runs_runtime_summary_table),by=num_code_sec)
diff_locations_runtime_indexes=seq(from=2,to=nrow(all_runs_runtime_summary_table),by=num_code_sec)
diff_areas_runtime_indexes=seq(from=3,to=nrow(all_runs_runtime_summary_table),by=num_code_sec)
diff_areas_alg_runtime_indexes=seq(from=4,to=nrow(all_runs_runtime_summary_table),by=num_code_sec)
motifs_one_diff_area_runtime_indexes=seq(from=5,to=nrow(all_runs_runtime_summary_table),by=num_code_sec)
seq_def_runtime_indexes=seq(from=6,to=nrow(all_runs_runtime_summary_table),by=num_code_sec)
pstree_modified_runtime_indexes=seq(from=7,to=nrow(all_runs_runtime_summary_table),by=num_code_sec)
get_nodes_to_keep_runtime_indexes=seq(from=8,to=nrow(all_runs_runtime_summary_table),by=num_code_sec)
# prune_updated_runtime_indexes=seq(from=9,to=nrow(all_runs_runtime_summary_table),by=num_code_sec)       if output plots!!!
motif_freq_pos_runtime_indexes=seq(from=9,to=nrow(all_runs_runtime_summary_table),by=num_code_sec)
predicting_all_diff_motifs_runtime_indexes=seq(from=10,to=nrow(all_runs_runtime_summary_table),by=num_code_sec)
removing_overlap_motifs_indexes=seq(from=11,to=nrow(all_runs_runtime_summary_table),by=num_code_sec)
diff_motifs_alg_runtime_indexes=seq(from=12,to=nrow(all_runs_runtime_summary_table),by=num_code_sec)
total_alg_runtime_indexes=seq(from=13,to=nrow(all_runs_runtime_summary_table),by=num_code_sec)

runtime_average_table_nrow=num_code_sec*num_measures

runtime_average_table=matrix(nrow=runtime_average_table_nrow)
first_col=c(rep("Average",num_code_sec),rep("Maximum",num_code_sec),rep("Minimum",num_code_sec),rep("Stand.Dev",num_code_sec))
runtime_average_table=cbind(runtime_average_table,first_col)

second_column=c(as.vector(all_runs_runtime_summary_table[1:num_code_sec,1]),as.vector(all_runs_runtime_summary_table[1:num_code_sec,1]),
                as.vector(all_runs_runtime_summary_table[1:num_code_sec,1]),as.vector(all_runs_runtime_summary_table[1:num_code_sec,1]))

runtime_average_table=cbind(runtime_average_table,second_column)

runtime_average_table=runtime_average_table[,-1]

row.names(runtime_average_table)=NULL

for(i in 2:ncol(all_runs_runtime_summary_table))
{
  current_column=all_runs_runtime_summary_table[,i]
  
  
  
#   average_sorting_data_runtime=round(mean(as.numeric(current_column[sorting_data_runtime_indexes])),digits=3)
#   average_generating_freq_table_runtime=round(mean(as.numeric(current_column[generating_freq_table_runtime_indexes])),digits=3)
#   average_generating_rel_freq_runtime=round(mean(as.numeric(current_column[generating_rel_freq_runtime_indexes])),digits=3)
  average_generating_pvalues_runtime=round(mean(as.numeric(current_column[generating_pvalues_runtime_indexes])),digits=3)
  average_diff_locations_runtime=round(mean(as.numeric(current_column[diff_locations_runtime_indexes])),digits=3)
  average_diff_areas_runtime=round(mean(as.numeric(current_column[diff_areas_runtime_indexes])),digits=3)
  
  average_diff_areas_alg_runtime=round(mean(as.numeric(current_column[diff_areas_alg_runtime_indexes])),digits=3)
  
  average_motifs_one_diff_area_runtime=round(mean(as.numeric(current_column[motifs_one_diff_area_runtime_indexes])),digits=3)
  average_seq_def_runtime=round(mean(as.numeric(current_column[seq_def_runtime_indexes])),digits=3)
  average_pstree_modified_runtime=round(mean(as.numeric(current_column[pstree_modified_runtime_indexes])),digits=3)
  average_get_nodes_to_keep_runtime=round(mean(as.numeric(current_column[get_nodes_to_keep_runtime_indexes])),digits=3)
#   average_prune_updated_runtime=round(mean(as.numeric(current_column[prune_updated_runtime_indexes])),digits=3)
  average_motif_freq_pos_runtime=round(mean(as.numeric(current_column[motif_freq_pos_runtime_indexes])),digits=3)
  
  average_predicting_all_diff_motifs_runtime=round(mean(as.numeric(current_column[predicting_all_diff_motifs_runtime_indexes])),digits=3)
  average_removing_overlap_motifs_runtime=round(mean(as.numeric(current_column[removing_overlap_motifs_indexes])),digits=3)
  
  average_diff_motifs_alg_runtime=round(mean(as.numeric(current_column[diff_motifs_alg_runtime_indexes])),digits=3)
  average_total_alg_runtime=round(mean(as.numeric(current_column[total_alg_runtime_indexes])),digits=3)
  
  
  
  
  
  
#   max_sorting_data_runtime=round(max(as.numeric(current_column[sorting_data_runtime_indexes])),digits=3)
#   max_generating_freq_table_runtime=round(max(as.numeric(current_column[generating_freq_table_runtime_indexes])),digits=3)
#   max_generating_rel_freq_runtime=round(max(as.numeric(current_column[generating_rel_freq_runtime_indexes])),digits=3)
  max_generating_pvalues_runtime=round(max(as.numeric(current_column[generating_pvalues_runtime_indexes])),digits=3)
  max_diff_locations_runtime=round(max(as.numeric(current_column[diff_locations_runtime_indexes])),digits=3)
  max_diff_areas_runtime=round(max(as.numeric(current_column[diff_areas_runtime_indexes])),digits=3)
  max_diff_areas_alg_runtime=round(max(as.numeric(current_column[diff_areas_alg_runtime_indexes])),digits=3)
  
  max_motifs_one_diff_area_runtime=round(max(as.numeric(current_column[motifs_one_diff_area_runtime_indexes])),digits=3)
  max_seq_def_runtime=round(max(as.numeric(current_column[seq_def_runtime_indexes])),digits=3)
  max_pstree_modified_runtime=round(max(as.numeric(current_column[pstree_modified_runtime_indexes])),digits=3)
  max_get_nodes_to_keep_runtime=round(max(as.numeric(current_column[get_nodes_to_keep_runtime_indexes])),digits=3)
#   max_prune_updated_runtime=round(max(as.numeric(current_column[prune_updated_runtime_indexes])),digits=3)
  max_motif_freq_pos_runtime=round(max(as.numeric(current_column[motif_freq_pos_runtime_indexes])),digits=3)
  
  max_predicting_all_diff_motifs_runtime=round(max(as.numeric(current_column[predicting_all_diff_motifs_runtime_indexes])),digits=3)
  max_removing_overlap_motifs_runtime=round(max(as.numeric(current_column[removing_overlap_motifs_indexes])),digits=3)
  
  max_diff_motifs_alg_runtime=round(max(as.numeric(current_column[diff_motifs_alg_runtime_indexes])),digits=3)
  max_total_alg_runtime=round(max(as.numeric(current_column[total_alg_runtime_indexes])),digits=3)
  
  
  
#   min_sorting_data_runtime=round(min(as.numeric(current_column[sorting_data_runtime_indexes])),digits=3)
#   min_generating_freq_table_runtime=round(min(as.numeric(current_column[generating_freq_table_runtime_indexes])),digits=3)
#   min_generating_rel_freq_runtime=round(min(as.numeric(current_column[generating_rel_freq_runtime_indexes])),digits=3)
  min_generating_pvalues_runtime=round(min(as.numeric(current_column[generating_pvalues_runtime_indexes])),digits=3)
  min_diff_locations_runtime=round(min(as.numeric(current_column[diff_locations_runtime_indexes])),digits=3)
  min_diff_areas_runtime=round(min(as.numeric(current_column[diff_areas_runtime_indexes])),digits=3)
  
  min_diff_areas_alg_runtime=round(min(as.numeric(current_column[diff_areas_alg_runtime_indexes])),digits=3)
  
  min_motifs_one_diff_area_runtime=round(min(as.numeric(current_column[motifs_one_diff_area_runtime_indexes])),digits=3)
  min_seq_def_runtime=round(min(as.numeric(current_column[seq_def_runtime_indexes])),digits=3)
  min_pstree_modified_runtime=round(min(as.numeric(current_column[pstree_modified_runtime_indexes])),digits=3)
  min_get_nodes_to_keep_runtime=round(min(as.numeric(current_column[get_nodes_to_keep_runtime_indexes])),digits=3)
#   min_prune_updated_runtime=round(min(as.numeric(current_column[prune_updated_runtime_indexes])),digits=3)
  min_motif_freq_pos_runtime=round(min(as.numeric(current_column[motif_freq_pos_runtime_indexes])),digits=3)
  
  
  min_predicting_all_diff_motifs_runtime=round(min(as.numeric(current_column[predicting_all_diff_motifs_runtime_indexes])),digits=3)
  min_removing_overlap_motifs_runtime=round(min(as.numeric(current_column[removing_overlap_motifs_indexes])),digits=3)
  
  
  min_diff_motifs_alg_runtime=round(min(as.numeric(current_column[diff_motifs_alg_runtime_indexes])),digits=3)
  min_total_alg_runtime=round(min(as.numeric(current_column[total_alg_runtime_indexes])),digits=3)
  
  
  
#   stdv_sorting_data_runtime=round(sd(as.numeric(current_column[sorting_data_runtime_indexes])),digits=3)
#   stdv_generating_freq_table_runtime=round(sd(as.numeric(current_column[generating_freq_table_runtime_indexes])),digits=3)
#   stdv_generating_rel_freq_runtime=round(sd(as.numeric(current_column[generating_rel_freq_runtime_indexes])),digits=3)
  stdv_generating_pvalues_runtime=round(sd(as.numeric(current_column[generating_pvalues_runtime_indexes])),digits=3)
  stdv_diff_locations_runtime=round(sd(as.numeric(current_column[diff_locations_runtime_indexes])),digits=3)
  stdv_diff_areas_runtime=round(sd(as.numeric(current_column[diff_areas_runtime_indexes])),digits=3)
  
  stdv_diff_areas_alg_runtime=round(sd(as.numeric(current_column[diff_areas_alg_runtime_indexes])),digits=3)
  
  stdv_motifs_one_diff_area_runtime=round(sd(as.numeric(current_column[motifs_one_diff_area_runtime_indexes])),digits=3)
  stdv_seq_def_runtime=round(sd(as.numeric(current_column[seq_def_runtime_indexes])),digits=3)
  stdv_pstree_modified_runtime=round(sd(as.numeric(current_column[pstree_modified_runtime_indexes])),digits=3)
  stdv_get_nodes_to_keep_runtime=round(sd(as.numeric(current_column[get_nodes_to_keep_runtime_indexes])),digits=3)
#   stdv_prune_updated_runtime=round(sd(as.numeric(current_column[prune_updated_runtime_indexes])),digits=3)
  stdv_motif_freq_pos_runtime=round(sd(as.numeric(current_column[motif_freq_pos_runtime_indexes])),digits=3)
  
  
  stdv_predicting_all_diff_motifs_runtime=round(sd(as.numeric(current_column[predicting_all_diff_motifs_runtime_indexes])),digits=3)
  stdv_removing_overlap_motifs_runtime=round(sd(as.numeric(current_column[removing_overlap_motifs_indexes])),digits=3)
  
  
  stdv_diff_motifs_alg_runtime=round(sd(as.numeric(current_column[diff_motifs_alg_runtime_indexes])),digits=3)
  stdv_total_alg_runtime=round(sd(as.numeric(current_column[total_alg_runtime_indexes])),digits=3)
  
  
  
  column_to_add=c(
                  average_generating_pvalues_runtime,
                  average_diff_locations_runtime,
                  average_diff_areas_runtime,
                  average_diff_areas_alg_runtime,
                  average_motifs_one_diff_area_runtime,
                  average_seq_def_runtime,
                  average_pstree_modified_runtime,
                  average_get_nodes_to_keep_runtime,
#                   average_prune_updated_runtime,
                  average_motif_freq_pos_runtime, 
                  average_predicting_all_diff_motifs_runtime,
                  average_removing_overlap_motifs_runtime,
                  average_diff_motifs_alg_runtime,
                  average_total_alg_runtime,
                  max_generating_pvalues_runtime,
                  max_diff_locations_runtime,
                  max_diff_areas_runtime,
                  max_diff_areas_alg_runtime,
                  max_motifs_one_diff_area_runtime,
                  max_seq_def_runtime,
                  max_pstree_modified_runtime,
                  max_get_nodes_to_keep_runtime,
#                   max_prune_updated_runtime,
                  max_motif_freq_pos_runtime,
                  max_predicting_all_diff_motifs_runtime,
                  max_removing_overlap_motifs_runtime,
                  max_diff_motifs_alg_runtime,
                  max_total_alg_runtime,
                  min_generating_pvalues_runtime,
                  min_diff_locations_runtime,
                  min_diff_areas_runtime,
                  min_diff_areas_alg_runtime,
                  min_motifs_one_diff_area_runtime,
                  min_seq_def_runtime,
                  min_pstree_modified_runtime,
                  min_get_nodes_to_keep_runtime,
#                   min_prune_updated_runtime,
                  min_motif_freq_pos_runtime,
                  min_predicting_all_diff_motifs_runtime,
                  min_removing_overlap_motifs_runtime,
                  min_diff_motifs_alg_runtime,
                  min_total_alg_runtime,
                  stdv_generating_pvalues_runtime,
                  stdv_diff_locations_runtime,
                  stdv_diff_areas_runtime,
                  stdv_diff_areas_alg_runtime,
                  stdv_motifs_one_diff_area_runtime,
                  stdv_seq_def_runtime,
                  stdv_pstree_modified_runtime,
                  stdv_get_nodes_to_keep_runtime,
#                   stdv_prune_updated_runtime,
                  stdv_motif_freq_pos_runtime,
                  stdv_predicting_all_diff_motifs_runtime,
                  stdv_removing_overlap_motifs_runtime,
                  stdv_diff_motifs_alg_runtime,
                  stdv_total_alg_runtime)
  
  runtime_average_table=cbind(runtime_average_table,column_to_add)
}



colnames(runtime_average_table)=c(" ",colnames(all_runs_runtime_summary_table))

runtime_avg_output_file=paste(data_type_dir,"/Table_runtime_summary_all_runs_MAX_MIN_AVG_STDV.csv",sep="")
write.csv(runtime_average_table,file=runtime_avg_output_file)

print("Changing format of table and plotting runtime as a function of variable")



variables=colnames(runtime_average_table)[-1:-2]
variables=unlist(strsplit(variables,"_",fixed=TRUE))
variable_name=variables[1]
index_to_keep=seq(from=0,to=length(variables),by=2)
variables=variables[index_to_keep]
num_variables=length(variables)
code_sections=second_column[1:num_code_sec]
# code_sections=c("Differentiable_Areas_Algorithm","Differentiable_Motifs_Algorithm","Complete_Algorithm")

max_col=c(as.numeric(runtime_average_table[length(code_sections)+(num_code_sec*1),3:(num_variables+2)]))
y_max_limit=max(max_col)

graph_table=data.frame()
for(i in 1:length(code_sections))
{
  code_section_col=c(rep(code_sections[i],times=num_variables))
  variable_col=variables
  variable_col=as.numeric(variable_col)
  
  avg_col=c(as.numeric(runtime_average_table[i,3:(num_variables+2)]))
  max_col=c(as.numeric(runtime_average_table[i+(num_code_sec*1),3:(num_variables+2)]))
  min_col=c(as.numeric(runtime_average_table[i+(num_code_sec*2),3:(num_variables+2)]))
  stdv_col=c(as.numeric(runtime_average_table[i+(num_code_sec*3),3:(num_variables+2)]))  
  graph_table=data.frame(Code_Section=factor(code_section_col),Variable_Parameter=factor(variable_col),
                         Average=as.numeric(avg_col),Maximum=as.numeric(max_col),Minimum=as.numeric(min_col),
                         Standard_Deviation=as.numeric(stdv_col))
  col_names=c("Code_Section","Variable_Parameter","Average","Maximum","Minimum","Standard_Deviation")
  colnames(graph_table)=col_names
  Variable_Parameter=reorder(graph_table$Variable_Parameter,variable_col)
  Code_Section=graph_table$Code_Section
  error_bar_min=graph_table$Average-graph_table$Standard_Deviation
  error_bar_min[which(error_bar_min<0)]=0
  error_bar_max=graph_table$Average+graph_table$Standard_Deviation
  #   graph_table_file=paste(data_type_dir,"/runtime_all_runs_MAX_MIN_AVG_STDV_Reformatted",i,".csv",sep="")
  #   write.csv(graph_table,file=graph_table_file)
  plots_dir=paste(data_type_dir,"/Graph_runtime_",i,"_",code_sections[i],".jpeg",sep="")
  plot=qplot(graph_table$Variable_Parameter,graph_table$Average,data=graph_table)
  plot=plot+geom_line(aes(group=Code_Section))
  plot=plot+geom_errorbar(aes(ymin=error_bar_min,ymax=error_bar_max,width=.1))
  plot=plot+ggtitle(paste(code_sections[i],"Runtime\nVS.",variable_name,"Over",num_runs,"Runs\n",param_list,data_type))
  plot=plot+theme(plot.title = element_text(size = 15))
  plot=plot+labs(x=variable_name,y="Average Runtime in Seconds")
  plot=plot+scale_x_discrete(limit = c(levels(graph_table$Variable_Parameter)))
  plot=plot+scale_y_continuous(limit = c(0,y_max_limit))
  #   plot=plot+scale_y_continuous(limit = c(min(graph_table$Minimum),max(graph_table$Maximum)))
  ggsave(plots_dir)
}
