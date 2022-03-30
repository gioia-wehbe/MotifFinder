



########################################################################################################
# Generating Average tables against all runs DIFFERENTIABLE MOTIFS PER LOCATION EVALUATION
########################################################################################################

param_list=getParamsForPlotTitle()


print("Generating evaluations on all runs for differentiable motifs per location")

#For differentiable motifs
pmin_values=getPminList()
num_pmin_values=length(pmin_values)

data_type_dir=""
if(data_type=="Training")
{
  data_type_dir=getTrainingDir()
}else{
  data_type_dir=getTestingDir()
}

all_motif_evaluations_dir=paste(data_type_dir,"/Table_predicted_motifs_per_location_metrics_all_runs.csv",sep="")
all_motifs_evaluations_table=read.csv(all_motif_evaluations_dir,header=FALSE,sep=",",,stringsAsFactors=FALSE)
all_motifs_evaluations_table=as.matrix(all_motifs_evaluations_table)
colnames(all_motifs_evaluations_table)=all_motifs_evaluations_table[1,]
all_motifs_evaluations_table=all_motifs_evaluations_table[-1,]

num_metrics=6

average_table_nrows=num_metrics*4*num_pmin_values    #num_metrics*num_av_max_min_stdv*num_pmin_values
motif_average_table=matrix(nrow=average_table_nrows)
nrow_per_stat=num_metrics*num_pmin_values
first_col=c(rep("Average",nrow_per_stat),rep("Maximum",nrow_per_stat),rep("Minimum",nrow_per_stat),rep("Stand.Dev",nrow_per_stat))
motif_average_table=cbind(motif_average_table,first_col)
second_2_columns=rbind(all_motifs_evaluations_table[1:nrow_per_stat,1:2],all_motifs_evaluations_table[1:nrow_per_stat,1:2],
                       all_motifs_evaluations_table[1:nrow_per_stat,1:2],all_motifs_evaluations_table[1:nrow_per_stat,1:2])
motif_average_table=cbind(motif_average_table,second_2_columns)
motif_average_table=motif_average_table[,-1]
row.names(motif_average_table)=NULL

for(c in 3:ncol(all_motifs_evaluations_table))    #for every column
{ 
  curr_column=all_motifs_evaluations_table[,c]
  
  
  average_values=c()
  max_values=c()
  min_values=c()
  std_dev_values=c()
  #   i=1
  start_ind=0
  for(i in 0:(num_pmin_values-1))   #for every pmin
  {
    start_ind=(i*num_metrics)+1
    current_pmin_sensitivity_indexes=seq(from=start_ind,to=nrow(all_motifs_evaluations_table),by=nrow_per_stat)
    current_pmin_specificity_indexes=seq(from=start_ind+1,to=nrow(all_motifs_evaluations_table),by=nrow_per_stat)
    current_pmin_precision_indexes=seq(from=start_ind+2,to=nrow(all_motifs_evaluations_table),by=nrow_per_stat)
    current_pmin_npv_indexes=seq(from=start_ind+3,to=nrow(all_motifs_evaluations_table),by=nrow_per_stat)
    current_pmin_ji_indexes=seq(from=start_ind+4,to=nrow(all_motifs_evaluations_table),by=nrow_per_stat)
    current_pmin_acc_indexes=seq(from=start_ind+5,to=nrow(all_motifs_evaluations_table),by=nrow_per_stat)
    
    current_pmin_avg_sen=mean(as.numeric(curr_column[current_pmin_sensitivity_indexes]))
    current_pmin_avg_spec=mean(as.numeric(curr_column[current_pmin_specificity_indexes]))
    current_pmin_avg_precision=mean(as.numeric(curr_column[current_pmin_precision_indexes]))
    current_pmin_avg_npv=mean(as.numeric(curr_column[current_pmin_npv_indexes]))
    current_pmin_avg_ji=mean(as.numeric(curr_column[current_pmin_ji_indexes]))
    current_pmin_avg_acc=mean(as.numeric(curr_column[current_pmin_acc_indexes]))
    average_values=c(average_values,current_pmin_avg_sen,current_pmin_avg_spec,current_pmin_avg_precision,
                     current_pmin_avg_npv,current_pmin_avg_ji,current_pmin_avg_acc)
    
    current_pmin_max_sen=as.numeric(max(curr_column[current_pmin_sensitivity_indexes]))
    current_pmin_max_spec=as.numeric(max(curr_column[current_pmin_specificity_indexes]))
    current_pmin_max_precision=as.numeric(max(curr_column[current_pmin_precision_indexes]))
    current_pmin_max_npv=as.numeric(max(curr_column[current_pmin_npv_indexes]))
    current_pmin_max_ji=as.numeric(max(curr_column[current_pmin_ji_indexes]))
    current_pmin_max_acc=as.numeric(max(curr_column[current_pmin_acc_indexes]))
    
    max_values=c(max_values,current_pmin_max_sen,current_pmin_max_spec,current_pmin_max_precision,
                 current_pmin_max_npv,current_pmin_max_ji,current_pmin_max_acc)
    
    current_pmin_min_sen=as.numeric(min(curr_column[current_pmin_sensitivity_indexes]))
    current_pmin_min_spec=as.numeric(min(curr_column[current_pmin_specificity_indexes]))
    current_pmin_min_precision=as.numeric(min(curr_column[current_pmin_precision_indexes]))
    current_pmin_min_npv=as.numeric(min(curr_column[current_pmin_npv_indexes]))
    current_pmin_min_ji=as.numeric(min(curr_column[current_pmin_ji_indexes]))
    current_pmin_min_acc=as.numeric(min(curr_column[current_pmin_acc_indexes]))
    min_values=c(min_values,current_pmin_min_sen,current_pmin_min_spec,current_pmin_min_precision,
                 current_pmin_min_npv,current_pmin_min_ji,current_pmin_min_acc)
    
    current_pmin_stdv_sen=as.numeric(sd(curr_column[current_pmin_sensitivity_indexes]))
    current_pmin_stdv_spec=as.numeric(sd(curr_column[current_pmin_specificity_indexes]))
    current_pmin_stdv_precision=as.numeric(sd(curr_column[current_pmin_precision_indexes]))
    current_pmin_stdv_npv=as.numeric(sd(curr_column[current_pmin_npv_indexes]))
    current_pmin_stdv_ji=as.numeric(sd(curr_column[current_pmin_ji_indexes]))
    current_pmin_stdv_acc=as.numeric(sd(curr_column[current_pmin_acc_indexes]))
    std_dev_values=c(std_dev_values,current_pmin_stdv_sen,current_pmin_stdv_spec,current_pmin_stdv_precision,
                     current_pmin_stdv_npv,current_pmin_stdv_ji,current_pmin_stdv_acc)
  } 
  average_curr_column=c(average_values,max_values,min_values,std_dev_values)
  motif_average_table=cbind(motif_average_table,average_curr_column)
}

colnames(motif_average_table)=c(" ",colnames(all_motifs_evaluations_table))

motif_avg_output_file=paste(data_type_dir,"/Table_predicted_motifs_per_locations_metrics_all_runs_MAX_MIN_AVG_STDV.csv",sep="")
write.csv(motif_average_table,file=motif_avg_output_file)



########################################################################################################
# Plotting Graphs FOR DIFFERENTIABLE MOTIFS PER LOCATION EVALUATION
########################################################################################################

print("Changing format of table and plotting for differentiable motifs per location evaluation")


motif_variables=colnames(motif_average_table)[-1:-3]
motif_variables=unlist(strsplit(motif_variables,"_",fixed=TRUE))
variable_name=motif_variables[1]
index_to_keep=seq(from=0,to=length(motif_variables),by=2)
motif_variables=motif_variables[index_to_keep]
num_motif_variables=length(motif_variables)
motif_variables=rep(motif_variables,times=num_pmin_values)

motif_graph_table=data.frame()
num_metrics=6

for( i in 1:num_metrics)
{
  
  pmin_col=c()
  for(pm in 1:num_pmin_values)
  {
    current_pmin_value=pmin_values[pm]
    pmin_repeated=rep(current_pmin_value,times=num_motif_variables)
    pmin_col=c(pmin_col,pmin_repeated)
  }
  variable_col=motif_variables
  variable_col=as.numeric(variable_col)
  metric=unlist(strsplit(as.character(motif_average_table[i,3]),"="))[1]
  
  avg_col=c()
  max_col=c()
  min_col=c()
  stdv_col=c()
  
  for(x in 0:(num_pmin_values-1))
  {
    ind=i+x*num_metrics
    avg_col=c(avg_col,motif_average_table[ind,4:(num_motif_variables+3)])
    ind=ind+nrow_per_stat
    max_col=c(max_col,motif_average_table[ind,4:(num_motif_variables+3)])
    ind=ind+nrow_per_stat
    min_col=c(min_col,motif_average_table[ind,4:(num_motif_variables+3)])
    ind=ind+nrow_per_stat
    stdv_col=c(stdv_col,motif_average_table[ind,4:(num_motif_variables+3)])
  }
  max_stdv=as.numeric(max(stdv_col))
  y_lim=1+max_stdv
  print("y_lim")
  print(y_lim)
  motif_graph_table=data.frame(Pmin_Value=factor(pmin_col),Variable_Parameter=factor(variable_col),Average=as.numeric(avg_col),
                               Maximum=as.numeric(max_col),Minimum=as.numeric(min_col),Standard_Deviation=as.numeric(stdv_col))
  col_names=c("Pmin_Value","Variable_Parameter","Average","Maximum","Minimum","Standard_Deviation")
  colnames(motif_graph_table)=col_names
  Variable_Parameter=reorder(motif_graph_table$Variable_Parameter,variable_col)
  #   average=motif_graph_table$Average
  #   average[which(average>1)]=1
  Pmin_Value=motif_graph_table$Pmin_Value
  error_bar_min=motif_graph_table$Average-motif_graph_table$Standard_Deviation
  error_bar_min[which(error_bar_min<0)]=0
  error_bar_max=motif_graph_table$Average+motif_graph_table$Standard_Deviation
  #   error_bar_max[which(error_bar_max>1)]=1
  plots_dir=paste(data_type_dir,"/Graph_location_evaluation_predicted_motifs_per_pmin_value_",metric,".jpeg",sep="")
  plot=qplot(motif_graph_table$Variable_Parameter,motif_graph_table$Average, colour=Pmin_Value, 
             data=motif_graph_table)
  plot=plot+geom_line(aes(group=motif_graph_table$Pmin_Value))
  plot=plot+geom_errorbar(aes(ymin=error_bar_min,ymax=error_bar_max,width=.1))
  plot=plot+ggtitle(paste(metric,"per location of Motifs\nVS.",variable_name,"Over",num_runs,"Runs\n",param_list,data_type)) 
  plot=plot+theme(plot.title = element_text(size =15))
  #   plot=plot+ggtitle(aes(label=param_list))
  plot=plot+labs(x=variable_name,y=paste("Average",metric))
  plot=plot+scale_x_discrete(limit = c(levels(motif_graph_table$Variable_Parameter)))
  plot=plot+scale_y_continuous(limit = c(0,y_lim))
  ggsave(plots_dir)
}

# motif_graph_table_file=paste(data_type_dir,"/predicted_motifs_per_location_per_pmin_value_metrics_all_runs_MAX_MIN_AVG_STDV_Reformatted.csv",sep="")
# write.csv(motif_graph_table,file=motif_graph_table_file)
