###############################################################
# Generating Average tables against all runs
###############################################################

OutputText=getOutputText()
OutputPlots=getOutputPlots()
param_list=getParamsForPlotTitle()

print("Generating SPLSDA evaluations on all runs...")


num_alleles=1 #ALL

num_metrics=6 #Spec, Sen, Acc, Prec, Recall,j-index
num_stats=4 #Avg Max Min Stdiv
num_runs=getNumRuns()

data_type_dir=""
if(data_type=="Training")
{
  data_type_dir=getTrainingAllOutputDir()
}else{
  data_type_dir=getTestingAllOutputDir()
}


#For differentiable locations
all_evaluations_dir=paste(data_type_dir,"/Table_splsda_differentiable_locations_per_allele_metrics_all_runs.csv",sep="")
all_evaluations_table=read.csv(all_evaluations_dir,header=FALSE,sep=",",stringsAsFactors=FALSE)

all_evaluations_table=as.matrix(all_evaluations_table)
colnames(all_evaluations_table)=all_evaluations_table[1,]
all_evaluations_table=all_evaluations_table[-1,]


average_table=matrix(nrow=num_alleles*num_metrics*num_stats)
print("nrow(average_table)")
print(nrow(average_table))
first_col=c(rep("Average",num_alleles*num_metrics),rep("Maximum",num_alleles*num_metrics),rep("Minimum",num_alleles*num_metrics),rep("Stand.Dev",num_alleles*num_metrics))
print("length(first_col)")
print(length(first_col))
average_table=cbind(average_table,first_col)
second_2_columns=c(all_evaluations_table[1:(num_alleles*num_metrics),3],all_evaluations_table[1:(num_alleles*num_metrics),3],
                       all_evaluations_table[1:(num_alleles*num_metrics),3],all_evaluations_table[1:(num_alleles*num_metrics),3])
print("nrow(second_2_columns)")
print(nrow(second_2_columns))
average_table=cbind(average_table,second_2_columns)
average_table=average_table[,-1]
row.names(average_table)=NULL


 
  all_SEN_indexes=seq(from=1,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
  all_SPC_indexes=seq(from=2,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
  all_PRC_indexes=seq(from=3,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
  all_NPV_indexes=seq(from=4,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
  all_JI_indexes=seq(from=5,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
  all_ACC_indexes=seq(from=6,to=nrow(all_evaluations_table),by=num_alleles*num_metrics)
  
  
  for(c in 4:ncol(all_evaluations_table))
  {
    curr_column=all_evaluations_table[,c]
    #   Computing Averages
    
    avg_all_SEN=mean(as.numeric(curr_column[all_SEN_indexes]),na.rm=TRUE)
    avg_all_SPC=mean(as.numeric(curr_column[all_SPC_indexes]),na.rm=TRUE)
    avg_all_PRC=mean(as.numeric(curr_column[all_PRC_indexes]),na.rm=TRUE)
    avg_all_NPV=mean(as.numeric(curr_column[all_NPV_indexes]),na.rm=TRUE)
    avg_all_JI=mean(as.numeric(curr_column[all_JI_indexes]),na.rm=TRUE)
    avg_all_ACC=mean(as.numeric(curr_column[all_ACC_indexes]),na.rm=TRUE)
    
    #   Computing Maximum
   
    max_all_SEN=as.numeric(max(curr_column[all_SEN_indexes],na.rm=TRUE))
    max_all_SPC=as.numeric(max(curr_column[all_SPC_indexes],na.rm=TRUE))
    max_all_PRC=as.numeric(max(curr_column[all_PRC_indexes],na.rm=TRUE))
    max_all_NPV=as.numeric(max(curr_column[all_NPV_indexes],na.rm=TRUE))
    max_all_JI=as.numeric(max(curr_column[all_JI_indexes],na.rm=TRUE))
    max_all_ACC=as.numeric(max(curr_column[all_ACC_indexes],na.rm=TRUE))
    
    #   Computing Minimum
    
    min_all_SEN=as.numeric(min(curr_column[all_SEN_indexes],na.rm=TRUE))
    min_all_SPC=as.numeric(min(curr_column[all_SPC_indexes],na.rm=TRUE))
    min_all_PRC=as.numeric(min(curr_column[all_PRC_indexes],na.rm=TRUE))
    min_all_NPV=as.numeric(min(curr_column[all_NPV_indexes],na.rm=TRUE))
    min_all_JI=as.numeric(min(curr_column[all_JI_indexes],na.rm=TRUE))
    min_all_ACC=as.numeric(min(curr_column[all_ACC_indexes],na.rm=TRUE))
    
    #   Computing Standard deviation
    
    sd_all_SEN=as.numeric(sd(curr_column[all_SEN_indexes],na.rm=TRUE))
    sd_all_SPC=as.numeric(sd(curr_column[all_SPC_indexes],na.rm=TRUE))
    sd_all_PRC=as.numeric(sd(curr_column[all_PRC_indexes],na.rm=TRUE))
    sd_all_NPV=as.numeric(sd(curr_column[all_NPV_indexes],na.rm=TRUE))
    sd_all_JI=as.numeric(sd(curr_column[all_JI_indexes],na.rm=TRUE))
    sd_all_ACC=as.numeric(sd(curr_column[all_ACC_indexes],na.rm=TRUE))
    
    average_curr_column=c(round(avg_all_SEN,digits=5),round(avg_all_SPC,digits=5),round(avg_all_PRC,digits=5)
                          ,round(avg_all_NPV,digits=5),round(avg_all_JI,digits=5),round(avg_all_ACC,digits=5)
                          ,round(max_all_SEN,digits=5),round(max_all_SPC,digits=5)
                          ,round(max_all_PRC,digits=5),round(max_all_NPV,digits=5),round(max_all_JI,digits=5)
                          ,round(max_all_ACC,digits=5),round(min_all_SEN,digits=5),round(min_all_SPC,digits=5)
                          ,round(min_all_PRC,digits=5),round(min_all_NPV,digits=5),round(min_all_JI,digits=5)
                          ,round(min_all_ACC,digits=5),round(sd_all_SEN,digits=5),round(sd_all_SPC,digits=5)
                          ,round(sd_all_PRC,digits=5),round(sd_all_NPV,digits=5),round(sd_all_JI,digits=5)
                          ,round(sd_all_ACC,digits=5))
    
    average_table=cbind(average_table,average_curr_column)
  }

colnames(average_table)=c(" ",colnames(all_evaluations_table)[3:ncol(all_evaluations_table)])

avg_output_file=paste(data_type_dir,"/Table_splsda_predicted_differentiable_locations_metrics_all_runs_MAX_MIN_AVG_STDV.csv",sep="")
write.csv(average_table,file=avg_output_file)



# #################################################################
# # Plotting Graphs for diff locations
# #################################################################
# 
# print("Changing format of table and plotting evaluation of differentiable locations prediction on all runs")
# 
# variables=colnames(average_table)[-1:-3]
# variables=variables[-length(variables)]
# variables=unlist(strsplit(variables,"_",fixed=TRUE))
# variable_name=variables[1]
# index_to_keep=seq(from=0,to=length(variables),by=2)
# variables=variables[index_to_keep]
# num_variables=length(variables)
# if(OutputText)
# {
#   variables=rep(variables,times=4)
# }
# 
# graph_table=data.frame()
# 
# num_metrics=6
# 
# if(OutputPlots)
# {
#   for( i in 1:num_metrics)
#   {
#     genotype_col=c(rep(0,times=num_variables),rep(1,times=num_variables),rep(2,times=num_variables),rep("all",times=num_variables))
#     variable_col=variables
#     variable_col=as.numeric(variable_col)
#     metric=unlist(strsplit(as.character(average_table[i,3]),"="))[1]
#     avg_col=as.numeric(c(average_table[i,4:(num_variables+3)],average_table[(i+5),4:(num_variables+3)],average_table[(i+10),4:(num_variables+3)],average_table[(i+15),4:(num_variables+3)]))
#     max_col=as.numeric(c(average_table[(i+20),4:(num_variables+3)],average_table[(i+25),4:(num_variables+3)],average_table[(i+30),4:(num_variables+3)],average_table[(i+35),4:(num_variables+3)]))
#     min_col=as.numeric(c(average_table[(i+40),4:(num_variables+3)],average_table[(i+45),4:(num_variables+3)],average_table[(i+50),4:(num_variables+3)],average_table[(i+55),4:(num_variables+3)]))
#     stdv_col=as.numeric(c(average_table[(i+60),4:(num_variables+3)],average_table[(i+65),4:(num_variables+3)],average_table[(i+70),4:(num_variables+3)],average_table[(i+75),4:(num_variables+3)]))
#     
#     graph_table=data.frame(Genotype=factor(genotype_col),Variable_Parameter=factor(variable_col),Average=avg_col,Maximum=max_col,Minimum=min_col,
#                            Standard_Deviation=stdv_col)
#     col_names=c("Genotype","Variable_Parameter","Average","Maximum","Minimum","Standard_Deviation")
#     colnames(graph_table)=col_names
#     Variable_Parameter=reorder(graph_table$Variable_Parameter,variable_col)
#     Genotype=graph_table$Genotype
#     error_bar_min=graph_table$Average-graph_table$Standard_Deviation
#     error_bar_min[which(error_bar_min<0)]=0
#     error_bar_max=graph_table$Average+graph_table$Standard_Deviation
#     error_bar_max[which(error_bar_max>1)]=1
#     plots_dir=paste(data_type_dir,"/Graph_predicted_differentiable_locations_per_allele_",metric,".jpeg",sep="")
#     plot=qplot(graph_table$Variable_Parameter,graph_table$Average, data=graph_table, colour=Genotype)
#     plot=plot+geom_line(aes(group=graph_table$Genotype))
#     plot=plot+geom_errorbar(aes(ymin=error_bar_min,ymax=error_bar_max,width=.1))
#     plot=plot+ggtitle(paste(metric,"of Differentiable SNP Locations\nVS.",variable_name,"Over",num_runs,"Runs\n",param_list))
#     plot=plot+theme(plot.title = element_text(size = 15))  
#     plot=plot+labs(x=variable_name,y=paste("Average",metric))
#     plot=plot+ scale_x_discrete(limit = c(levels(graph_table$Variable_Parameter)))
#     plot=plot+scale_y_continuous(limit=c(0,1))
#     ggsave(plots_dir)
#   }
# }else {
#   for( i in 1:num_metrics)
#   {
#     genotype_col=c(rep("all",times=num_variables))
#     print("num_variables")
#     print(num_variables)
#     print("variables")
#     print(variables)
#     variable_col=variables
#     variable_col=as.numeric(variable_col)
#     metric=unlist(strsplit(as.character(average_table[i,3]),"="))[1]
#     avg_col=as.numeric(c(average_table[i,4:(num_variables+3)]))
#     max_col=as.numeric(c(average_table[(i+6),4:(num_variables+3)]))
#     min_col=as.numeric(c(average_table[(i+12),4:(num_variables+3)]))
#     stdv_col=as.numeric(c(average_table[(i+18),4:(num_variables+3)]))
#     
#     graph_table=data.frame(Genotype=factor(genotype_col),Variable_Parameter=factor(variable_col),Average=avg_col,Maximum=max_col,Minimum=min_col,
#                            Standard_Deviation=stdv_col)
#     col_names=c("Genotype","Variable_Parameter","Average","Maximum","Minimum","Standard_Deviation")
#     colnames(graph_table)=col_names
# #     graph_table_file=paste(data_type_dir,"/predicted_differentiable_locations_per_allele_metrics_all_runs_MAX_MIN_AVG_STDV_Reformatted",i,".csv",sep="")
# #     write.csv(graph_table,file=graph_table_file)
#     Variable_Parameter=reorder(graph_table$Variable_Parameter,variable_col)
#     Genotype=graph_table$Genotype
#     error_bar_min=graph_table$Average-graph_table$Standard_Deviation
#     error_bar_min[which(error_bar_min<0)]=0
#     error_bar_max=graph_table$Average+graph_table$Standard_Deviation
#     error_bar_max[which(error_bar_max>1)]=1
#     plots_dir=paste(data_type_dir,"/Graph_predicted_differentiable_locations_per_allele_",metric,".jpeg",sep="")
#     plot=qplot(graph_table$Variable_Parameter,graph_table$Average, data=graph_table, colour=Genotype)
#     plot=plot+geom_line(aes(group=graph_table$Genotype))
#     plot=plot+geom_errorbar(aes(ymin=error_bar_min,ymax=error_bar_max,width=.1))
#     plot=plot+ggtitle(paste(metric,"of Differentiable SNP Locations\nVS.",variable_name,"Over",num_runs,"Runs\n",param_list))
#     plot=plot+theme(plot.title = element_text(size = 15))  
#     plot=plot+labs(x=variable_name,y=paste("Average",metric))
#     plot=plot+ scale_x_discrete(limit = c(levels(graph_table$Variable_Parameter)))
#     plot=plot+scale_y_continuous(limit=c(0,1))
#     ggsave(plots_dir)
#   }
# }
# 
# 
# 
# 
# 
# 
