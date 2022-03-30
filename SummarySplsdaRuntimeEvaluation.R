#########################################
# Plotting runtime graph
#########################################


print("Plotting SPLSDA Runtime Graph")

param_list=getParamsForPlotTitle()

#read table of all runs

data_type_dir=""
if(data_type=="Training")
{
  data_type_dir=getTrainingAllOutputDir()
}else if(data_type=="Testing"){
  data_type_dir=getTestingAllOutputDir()
}else {
  data_type_dir=getCompleteDataRuntimeEvaluationDir()
}

all_runs_runtime_summary_file=paste(data_type_dir,"/Table_splsda_runtime_summary_all_runs.csv",sep="")
all_runs_runtime_summary_table=read.csv(all_runs_runtime_summary_file,sep=",",header=FALSE, colClasses="character",stringsAsFactors=FALSE)
col_names=all_runs_runtime_summary_table[1,]
colnames(all_runs_runtime_summary_table)=col_names
all_runs_runtime_summary_table=all_runs_runtime_summary_table[-1,]

num_code_sec=3
num_measures=4 #Avg,max,min,stdv

running_splsda_runtime_indexes=seq(from=1,to=nrow(all_runs_runtime_summary_table),by=num_code_sec)
selecting_vars_splsda_runtime_indexes=seq(from=2,to=nrow(all_runs_runtime_summary_table),by=num_code_sec)
total_splsda_runtime_indexes=seq(from=3,to=nrow(all_runs_runtime_summary_table),by=num_code_sec)

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
  
  average_running_splsda_runtime=round(mean(as.numeric(current_column[running_splsda_runtime_indexes])),digits=3)
  average_selecting_vars_splsda_runtime=round(mean(as.numeric(current_column[selecting_vars_splsda_runtime_indexes])),digits=3)
  average_total_splsda_runtime=round(mean(as.numeric(current_column[total_splsda_runtime_indexes])),digits=3)
  
  max_running_splsda_runtime=round(max(as.numeric(current_column[running_splsda_runtime_indexes])),digits=3)
  max_selecting_vars_splsda_runtime=round(max(as.numeric(current_column[selecting_vars_splsda_runtime_indexes])),digits=3)
  max_total_splsda_runtime=round(max(as.numeric(current_column[total_splsda_runtime_indexes])),digits=3)
  
  min_running_splsda_runtime=round(min(as.numeric(current_column[running_splsda_runtime_indexes])),digits=3)
  min_selecting_vars_splsda_runtime=round(min(as.numeric(current_column[selecting_vars_splsda_runtime_indexes])),digits=3)
  min_total_splsda_runtime=round(min(as.numeric(current_column[total_splsda_runtime_indexes])),digits=3)
  
  stdv_running_splsda_runtime=round(sd(as.numeric(current_column[running_splsda_runtime_indexes])),digits=3)
  stdv_selecting_vars_splsda_runtime=round(sd(as.numeric(current_column[selecting_vars_splsda_runtime_indexes])),digits=3)
  stdv_total_splsda_runtime=round(sd(as.numeric(current_column[total_splsda_runtime_indexes])),digits=3)
  
  column_to_add=c(average_running_splsda_runtime,average_selecting_vars_splsda_runtime,average_total_splsda_runtime,
                  max_running_splsda_runtime,max_selecting_vars_splsda_runtime,max_total_splsda_runtime,
                  min_running_splsda_runtime,min_selecting_vars_splsda_runtime,min_total_splsda_runtime,
                  stdv_running_splsda_runtime,stdv_selecting_vars_splsda_runtime,stdv_total_splsda_runtime)
  
  runtime_average_table=cbind(runtime_average_table,column_to_add)
}



colnames(runtime_average_table)=c(" ",colnames(all_runs_runtime_summary_table))

runtime_avg_output_file=paste(data_type_dir,"/Table_splsda_runtime_summary_all_runs_MAX_MIN_AVG_STDV.csv",sep="")
write.csv(runtime_average_table,file=runtime_avg_output_file)



print("Changing format of table and plotting splsda runtime as a function of variables")


alg_all_runs_runtime_summary_file=paste(data_type_dir,"/Table_runtime_summary_all_runs_MAX_MIN_AVG_STDV.csv",sep="")
alg_runtime_average_table=read.csv(alg_all_runs_runtime_summary_file,sep=",",header=FALSE, colClasses="character",stringsAsFactors=FALSE)
col_names=alg_runtime_average_table[1,]
colnames(alg_runtime_average_table)=col_names
alg_runtime_average_table=alg_runtime_average_table[-1,]
alg_num_code_sec=14


print("head(alg_runtime_average_table)")
print(head(alg_runtime_average_table))

variables=colnames(runtime_average_table)[-1:-2]
variables=unlist(strsplit(variables,"_",fixed=TRUE))
variable_name=variables[1]
index_to_keep=seq(from=0,to=length(variables),by=2)
variables=variables[index_to_keep]
num_variables=length(variables)
code_sections=second_column[1:num_code_sec]
# code_sections=c("Differentiable_Areas_Algorithm","Differentiable_Motifs_Algorithm","Complete_Algorithm")

max_col=c(as.numeric(runtime_average_table[length(code_sections)+(num_code_sec*1),3:(num_variables+2)]),
          as.numeric(alg_runtime_average_table[4+(alg_num_code_sec*1),4:(num_variables+3)]))
print("max_col")
print(max_col)
y_max_limit=as.numeric(max(max_col))
print("y_max_limit")
print(y_max_limit)

  graph_table=data.frame()
  alg_col=c(rep("SPLSDA",num_variables),rep("SNP-DMF",num_variables))
  code_section_col=c(rep("Feature Selection Runtime Benchmarks",times=num_variables*2))
  variable_col=variables
  variable_col=as.numeric(variable_col)
  variable_col=rep(variable_col,2)
  avg_col=c(as.numeric(runtime_average_table[3,3:(num_variables+2)]),as.numeric(alg_runtime_average_table[4,4:(num_variables+3)]))
  max_col=c(as.numeric(runtime_average_table[3+(num_code_sec*1),3:(num_variables+2)]),as.numeric(alg_runtime_average_table[4+(alg_num_code_sec*1),4:(num_variables+3)]))
  min_col=c(as.numeric(runtime_average_table[3+(num_code_sec*2),3:(num_variables+2)]),as.numeric(alg_runtime_average_table[4+(alg_num_code_sec*2),4:(num_variables+3)]))
  stdv_col=c(as.numeric(runtime_average_table[3+(num_code_sec*3),3:(num_variables+2)]),as.numeric(alg_runtime_average_table[4+(alg_num_code_sec*3),4:(num_variables+3)]))  
max_stdv=as.numeric(max(stdv_col)) 
graph_table=data.frame(Algorithm=factor(alg_col),Variable_Parameter=factor(variable_col),
                         Average=as.numeric(avg_col),Maximum=as.numeric(max_col),Minimum=as.numeric(min_col),
                         Standard_Deviation=as.numeric(stdv_col))
  col_names=c("Algorithm","Variable_Parameter","Average","Maximum","Minimum","Standard_Deviation")
  colnames(graph_table)=col_names
  Variable_Parameter=reorder(graph_table$Variable_Parameter,variable_col)
  Algorithm=graph_table$Algorithm
  error_bar_min=graph_table$Average-graph_table$Standard_Deviation
  error_bar_min[which(error_bar_min<0)]=0
  error_bar_max=graph_table$Average+graph_table$Standard_Deviation
#   graph_table_file=paste(data_type_dir,"/Table_runtime_all_runs_benchmarks_MAX_MIN_AVG_STDV_Reformatted.csv",sep="")
#   write.csv(graph_table,file=graph_table_file)
  plots_dir=paste(data_type_dir,"/Graph_runtime_",code_section_col[1],".jpeg",sep="")
  plot=qplot(graph_table$Variable_Parameter,graph_table$Average,data=graph_table,colour=Algorithm)
  plot=plot+geom_line(aes(group=Algorithm))
  plot=plot+geom_errorbar(aes(ymin=error_bar_min,ymax=error_bar_max,width=.1))
  plot=plot+ggtitle(paste("Benchmarking Runtime Over",num_runs,"Runs\n",param_list,data_type))
  plot=plot+theme(plot.title = element_text(size = 15))
  plot=plot+labs(x=variable_name,y="Average Runtime in Seconds")
  plot=plot+scale_x_discrete(limit = c(levels(graph_table$Variable_Parameter)))
  y_max_limit=y_max_limit+max_stdv
  plot=plot+scale_y_continuous(limit = c(0,y_max_limit))
  #   plot=plot+scale_y_continuous(limit = c(min(graph_table$Minimum),max(graph_table$Maximum)))
  ggsave(plots_dir)
