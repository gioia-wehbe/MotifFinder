
##########################################################################
# Generating Average tables against all runs DIFFERENTIABLE MOTIFS
##########################################################################

print("inside SummaryFullDiffMotifsEvaluation.R")

param_list=getParamsForPlotTitle()
num_folds=getNumFolds()

print("Generating evaluations on all runs for differentiable motifs")

num_metrics=6

for(dt in 1:length(data_type_vector))
{
  data_type=data_type_vector[dt]
  print("data_type")
  print(data_type)
  data_type_dir=""
  if(data_type=="Training")
  {
    data_type_dir=getTrainingAllOutputDir()
  }else{
    data_type_dir=getTestingAllOutputDir()
  }
  
  
  #For differentiable motifs
  pmin_values=getPminList()
  num_pmin_values=length(pmin_values)
  
  print("data_type_dir")
  print(data_type_dir)
  
  
  all_motif_evaluations_dir=paste(data_type_dir,"/Table_predicted_motifs_complete_metrics_all_runs.csv",sep="")
  all_motifs_evaluations_table=read.csv(all_motif_evaluations_dir,header=FALSE,sep=",",,stringsAsFactors=FALSE)
  all_motifs_evaluations_table=as.matrix(all_motifs_evaluations_table)
  colnames(all_motifs_evaluations_table)=all_motifs_evaluations_table[1,]
  all_motifs_evaluations_table=all_motifs_evaluations_table[-1,]
  
  print("all_motif_evaluations_dir")
  print(all_motif_evaluations_dir)
  print("head(all_motifs_evaluations_table)")
  print(head(all_motifs_evaluations_table))
  
  num_metrics=6
  average_table_nrows=num_metrics*4*num_pmin_values    #num_metrics*num_av_max_min_stdv*num_pmin_values
  motif_average_table=matrix(nrow=average_table_nrows)
  nrow_per_stat=num_metrics*num_pmin_values
  first_col=c(rep("Average",nrow_per_stat),rep("Maximum",nrow_per_stat),rep("Minimum",nrow_per_stat),rep("Stand.Dev",nrow_per_stat))
  motif_average_table=cbind(motif_average_table,first_col)
  second_2_columns=rbind(all_motifs_evaluations_table[1:nrow_per_stat,3:4],all_motifs_evaluations_table[1:nrow_per_stat,3:4],
                         all_motifs_evaluations_table[1:nrow_per_stat,3:4],all_motifs_evaluations_table[1:nrow_per_stat,3:4])
  motif_average_table=cbind(motif_average_table,second_2_columns)
  motif_average_table=motif_average_table[,-1]
  row.names(motif_average_table)=NULL
  
  for(c in 5:ncol(all_motifs_evaluations_table))    #for every column
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
      
      current_pmin_SEN_indexes=seq(from=start_ind,to=nrow(all_motifs_evaluations_table),by=nrow_per_stat)
      current_pmin_SPC_indexes=seq(from=start_ind+1,to=nrow(all_motifs_evaluations_table),by=nrow_per_stat)
      current_pmin_PRC_indexes=seq(from=start_ind+2,to=nrow(all_motifs_evaluations_table),by=nrow_per_stat)
      current_pmin_NPV_indexes=seq(from=start_ind+3,to=nrow(all_motifs_evaluations_table),by=nrow_per_stat)
      current_pmin_JI_indexes=seq(from=start_ind+4,to=nrow(all_motifs_evaluations_table),by=nrow_per_stat)
      current_pmin_ACC_indexes=seq(from=start_ind+5,to=nrow(all_motifs_evaluations_table),by=nrow_per_stat)
      
      current_pmin_avg_SEN=mean(as.numeric(curr_column[current_pmin_SEN_indexes]),na.rm=TRUE)
      current_pmin_avg_SPC=mean(as.numeric(curr_column[current_pmin_SPC_indexes]),na.rm=TRUE)
      current_pmin_avg_PRC=mean(as.numeric(curr_column[current_pmin_PRC_indexes]),na.rm=TRUE)
      current_pmin_avg_NPV=mean(as.numeric(curr_column[current_pmin_NPV_indexes]),na.rm=TRUE)
      current_pmin_avg_JI=mean(as.numeric(curr_column[current_pmin_JI_indexes]),na.rm=TRUE)
      current_pmin_avg_ACC=mean(as.numeric(curr_column[current_pmin_ACC_indexes]),na.rm=TRUE)
      average_values=c(average_values,current_pmin_avg_SEN,current_pmin_avg_SPC,current_pmin_avg_PRC,current_pmin_avg_NPV,
                       current_pmin_avg_JI,current_pmin_avg_ACC)
      
      current_pmin_max_SEN=as.numeric(max(curr_column[current_pmin_SEN_indexes],na.rm=TRUE))
      current_pmin_max_SPC=as.numeric(max(curr_column[current_pmin_SPC_indexes],na.rm=TRUE))
      current_pmin_max_PRC=as.numeric(max(curr_column[current_pmin_PRC_indexes],na.rm=TRUE))
      current_pmin_max_NPV=as.numeric(max(curr_column[current_pmin_NPV_indexes],na.rm=TRUE))
      current_pmin_max_JI=as.numeric(max(curr_column[current_pmin_JI_indexes],na.rm=TRUE))
      current_pmin_max_ACC=as.numeric(max(curr_column[current_pmin_ACC_indexes],na.rm=TRUE))
      max_values=c(max_values,current_pmin_max_SEN,current_pmin_max_SPC,current_pmin_max_PRC,current_pmin_max_NPV,
                   current_pmin_max_JI,current_pmin_max_ACC)
      
      current_pmin_min_SEN=as.numeric(min(curr_column[current_pmin_SEN_indexes],na.rm=TRUE))
      current_pmin_min_SPC=as.numeric(min(curr_column[current_pmin_SPC_indexes],na.rm=TRUE))
      current_pmin_min_PRC=as.numeric(min(curr_column[current_pmin_PRC_indexes],na.rm=TRUE))
      current_pmin_min_NPV=as.numeric(min(curr_column[current_pmin_NPV_indexes],na.rm=TRUE))
      current_pmin_min_JI=as.numeric(min(curr_column[current_pmin_JI_indexes],na.rm=TRUE))
      current_pmin_min_ACC=as.numeric(min(curr_column[current_pmin_ACC_indexes],na.rm=TRUE))
      min_values=c(min_values,current_pmin_min_SEN,current_pmin_min_SPC,current_pmin_min_PRC,current_pmin_min_NPV
                   ,current_pmin_min_JI,current_pmin_min_ACC)
      
      current_pmin_stdv_SEN=as.numeric(sd(curr_column[current_pmin_SEN_indexes],na.rm=TRUE))
      current_pmin_stdv_SPC=as.numeric(sd(curr_column[current_pmin_SPC_indexes],na.rm=TRUE))
      current_pmin_stdv_PRC=as.numeric(sd(curr_column[current_pmin_PRC_indexes],na.rm=TRUE))
      current_pmin_stdv_NPV=as.numeric(sd(curr_column[current_pmin_NPV_indexes],na.rm=TRUE))
      current_pmin_stdv_JI=as.numeric(sd(curr_column[current_pmin_JI_indexes],na.rm=TRUE))
      current_pmin_stdv_ACC=as.numeric(sd(curr_column[current_pmin_ACC_indexes],na.rm=TRUE))

      
      
      
      std_dev_values=c(std_dev_values,current_pmin_stdv_SEN,current_pmin_stdv_SPC,current_pmin_stdv_PRC,current_pmin_stdv_NPV
                       ,current_pmin_stdv_JI,current_pmin_stdv_ACC)
      
    } 
    average_curr_column=c(average_values,max_values,min_values,std_dev_values)
    motif_average_table=cbind(motif_average_table,average_curr_column)
  }
  
  
  all_motifs_evaluations_table_colnames=colnames(all_motifs_evaluations_table)
  colnames(motif_average_table)=c(" ",all_motifs_evaluations_table_colnames[3:length(all_motifs_evaluations_table_colnames)])
  
  motif_avg_output_file=paste(data_type_dir,"/Table_predicted_motifs_complete_metrics_all_runs_MAX_MIN_AVG_STDV.csv",sep="")
  write.csv(motif_average_table,file=motif_avg_output_file)
  print("motif_avg_output_file")
  print(motif_avg_output_file)
  
  
}#end of for(dt in 1:length(data_type_vector))






#################################################################
# Plotting Graphs FOR DIFFERENTIABLE MOTIFS
#################################################################

print("Changing format of table and plotting for differentiable motifs")



#GET TRAINING TABLES

data_type="Training"
data_type_dir=getTrainingAllOutputDir()

average_table_train_file=paste(data_type_dir,"/Table_predicted_motifs_complete_metrics_all_runs_MAX_MIN_AVG_STDV.csv",sep="")
average_table_train=read.csv(average_table_train_file,header=FALSE,sep=",",stringsAsFactors=FALSE)
print("average_table_train_file")
print(average_table_train_file)


average_table_train=as.matrix(average_table_train)
colnames(average_table_train)=average_table_train[1,]
average_table_train=average_table_train[-1,]
average_table_train=average_table_train[,-1]

print("head(average_table_train)")
print(head(average_table_train))




#GET TESTING TABLES

data_type="Testing"
data_type_dir=getTestingAllOutputDir()

average_table_test_file=paste(data_type_dir,"/Table_predicted_motifs_complete_metrics_all_runs_MAX_MIN_AVG_STDV.csv",sep="")
average_table_test=read.csv(average_table_test_file,header=FALSE,sep=",",stringsAsFactors=FALSE)

print("average_table_test_file")
print(average_table_test_file)

average_table_test=as.matrix(average_table_test)
colnames(average_table_test)=average_table_test[1,]
average_table_test=average_table_test[-1,]
average_table_test=average_table_test[,-1]

print("head(average_table_test)")
print(head(average_table_test))










motif_variables=colnames(average_table_train)[-1:-3]
motif_variables=unlist(strsplit(motif_variables,"_",fixed=TRUE))
variable_name=motif_variables[1]
index_to_keep=seq(from=0,to=length(motif_variables),by=2)
motif_variables=motif_variables[index_to_keep]
num_motif_variables=length(motif_variables)
motif_variables=rep(motif_variables,times=num_pmin_values)

motif_graph_table=data.frame()
# num_metrics=3

all_stdv_col=as.numeric(c(average_table_train[(1+18),4:(num_motif_variables+3)],
                          average_table_train[(2+18),4:(num_motif_variables+3)],
                          average_table_train[(3+18),4:(num_motif_variables+3)],
                          average_table_train[(4+18),4:(num_motif_variables+3)],
                          average_table_train[(5+18),4:(num_motif_variables+3)],
                          average_table_train[(6+18),4:(num_motif_variables+3)],
                          average_table_test[(1+18),4:(num_motif_variables+3)],
                          average_table_test[(2+18),4:(num_motif_variables+3)],
                          average_table_test[(3+18),4:(num_motif_variables+3)],
                          average_table_test[(4+18),4:(num_motif_variables+3)],
                          average_table_test[(5+18),4:(num_motif_variables+3)],
                          average_table_test[(6+18),4:(num_motif_variables+3)]))



print("all_stdv_col")
print(all_stdv_col)
all_stdv_col
max_stdv=as.numeric(max(all_stdv_col,na.rm=TRUE))
print("max_stdv")
print(max_stdv)
y_lim=1+max_stdv
print("y_lim")
print(y_lim)

for( i in 1:num_metrics)
{
  alg_col=rep("SNP-DMF",num_motif_variables*2)
  data_type_col=c(rep("Training",num_motif_variables),rep("Testing",num_motif_variables))
  print("data_type_col")
  print(data_type_col)
  
  pmin_col=c()
  for(pm in 1:num_pmin_values)
  {
    current_pmin_value=pmin_values[pm]
    pmin_repeated=rep(current_pmin_value,times=num_motif_variables)
    pmin_col=c(pmin_col,pmin_repeated)
  }
  pmin_col=rep(pmin_col,2)
  print("pmin_col")
  print(pmin_col)
  
  variable_col=motif_variables
  variable_col=as.numeric(variable_col)
  variable_col=rep(variable_col,2)
  print("variable_col")
  print(variable_col)
  
  metric=unlist(strsplit(as.character(motif_average_table[i,3]),"="))[1]
  
  avg_col=c()
  max_col=c()
  min_col=c()
  stdv_col=c()
  
  for(x in 0:(num_pmin_values-1))
  {
    ind=i+x*num_metrics
    avg_col=c(avg_col,average_table_train[ind,4:(num_motif_variables+3)],average_table_test[ind,4:(num_motif_variables+3)])
    ind=ind+nrow_per_stat
    max_col=c(max_col,average_table_train[ind,4:(num_motif_variables+3)],average_table_test[ind,4:(num_motif_variables+3)])
    ind=ind+nrow_per_stat
    min_col=c(min_col,average_table_train[ind,4:(num_motif_variables+3)],average_table_test[ind,4:(num_motif_variables+3)])
    ind=ind+nrow_per_stat
    stdv_col=c(stdv_col,average_table_train[ind,4:(num_motif_variables+3)],average_table_test[ind,4:(num_motif_variables+3)])
  }
  
  print("avg_col")
  print(avg_col)
  print("max_col")
  print(max_col)
  print("min_col")
  print(min_col)
  print("stdv_col")
  print(stdv_col)
  
  
  motif_graph_table=data.frame(Algorithm=factor(alg_col),Data_Type=factor(data_type_col),Pmin_Value=factor(pmin_col),Variable_Parameter=factor(variable_col),Average=as.numeric(avg_col),
                               Maximum=as.numeric(max_col),Minimum=as.numeric(min_col),Standard_Deviation=as.numeric(stdv_col))
  col_names=c("Algorithm","Data_Type","Pmin_Value","Variable_Parameter","Average","Maximum","Minimum","Standard_Deviation")
  colnames(motif_graph_table)=col_names
#   graph_table_file=paste(getTrainingTestingPerformanceEvaluationDir(),"/Table_predicted_motifs_complete_metrics_all_runs_MAX_MIN_AVG_STDV_Reformatted",metric,".csv",sep="")
#   write.csv(motif_graph_table,file=graph_table_file)
  
  
  #################################################################
  # Plotting Graphs FOR DIFFERENTIABLE MOTIFS
  #################################################################
  
  print("Plotting...")
  Algorithm=motif_graph_table$Algorithm
  Data_Type=motif_graph_table$Data_Type
  Variable_Parameter=reorder(motif_graph_table$Variable_Parameter,variable_col)
  #   average=motif_graph_table$Average
  #   average[which(average>1)]=1
  Pmin_Value=motif_graph_table$Pmin_Value
  error_bar_min=motif_graph_table$Average-motif_graph_table$Standard_Deviation
  error_bar_min[which(error_bar_min<0)]=0
  error_bar_max=motif_graph_table$Average+motif_graph_table$Standard_Deviation
  plots_dir=paste(getTrainingTestingPerformanceEvaluationDir(),"/Graph_complete_evaluation_predicted_motifs_per_pmin_value_",metric,".jpeg",sep="")
#   plot=qplot(motif_graph_table$Variable_Parameter,motif_graph_table$Average, colour=Pmin_Value,data=motif_graph_table)
  plot=qplot(motif_graph_table$Variable_Parameter,motif_graph_table$Average,data=motif_graph_table,colour=Algorithm)
  plot=plot+geom_line(aes(group=motif_graph_table$Data_Type,linetype=Data_Type))
  plot=plot+geom_errorbar(aes(ymin=error_bar_min,ymax=error_bar_max,width=.1))
  plot=plot+ggtitle(paste(metric,"of Motifs\nVS.",variable_name,"Over",num_runs,"Runs","and ",num_folds,"Folds\n",param_list)) 
  plot=plot+theme(plot.title = element_text(size =15))
  #   plot=plot+ggtitle(aes(label=param_list))
  plot=plot+labs(x=variable_name,y=paste("Average",metric))
  plot=plot+scale_x_discrete(limit = c(levels(motif_graph_table$Variable_Parameter)))
  plot=plot+scale_y_continuous(limit = c(0,y_lim))
  ggsave(plots_dir)
}

# motif_graph_table_file=paste(data_type_dir,"/predicted_motifs_per_pmin_value_metrics_all_runs_MAX_MIN_AVG_STDV_Reformatted.csv",sep="")
# write.csv(motif_graph_table,file=motif_graph_table_file)

